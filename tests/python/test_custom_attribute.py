"""Tests for SPZ_NIANTIC_custom_attribute and the zstd skippable frames that carry its values."""

import struct

import numpy as np
import pytest

import spz

if not spz.has_extension_support():
    pytest.skip("built without extension support", allow_module_level=True)

S = spz.SpzAttributeSemantic
HEADER = struct.Struct("<IIIBBBBI12x")
SKIPPABLE_MAGIC = 0x184D2A53
ZSTD_MAGIC = 0xFD2FB528
CUSTOM_ATTRIBUTE_TYPE = 0x4E530002
FLAG_HAS_EXTENSIONS = 0x2
MAX_RAW_BLOCK = 128 * 1024
# Semantic values as stored in the record.
CUSTOM, SPATIAL = 0, 1


def _make_cloud(num_points=50, sh_degree=1, seed=0):
    rng = np.random.default_rng(seed)
    cloud = spz.GaussianCloud()
    cloud.sh_degree = sh_degree
    cloud.positions = rng.uniform(-1.0, 1.0, size=num_points * 3).astype(np.float32)
    cloud.scales = rng.uniform(-2.0, 2.0, size=num_points * 3).astype(np.float32)
    cloud.rotations = rng.uniform(-1.0, 1.0, size=num_points * 4).astype(np.float32)
    cloud.alphas = rng.uniform(-1.0, 1.0, size=num_points).astype(np.float32)
    cloud.colors = rng.uniform(0.0, 1.0, size=num_points * 3).astype(np.float32)
    sh_dim = {0: 0, 1: 3, 2: 8, 3: 15}[sh_degree]
    cloud.sh = rng.uniform(-0.5, 0.5, size=num_points * sh_dim * 3).astype(np.float32)
    return cloud


def _attribute(name, dim, values, bits=32, semantic=S.CUSTOM):
    attr = spz.SpzExtensionCustomAttributeNiantic()
    attr.name = name
    attr.dim = dim
    attr.semantic = semantic
    attr.bits = bits
    attr.values = np.asarray(values, dtype=np.float32)
    return attr


def _decode_all():
    options = spz.UnpackOptions()
    options.custom_attributes = None
    return options


def _round_trip(cloud, pack=None, unpack=None):
    data = spz.save_spz_to_buffer(cloud, pack or spz.PackOptions())
    return data, spz.load_spz_from_buffer(data, unpack or _decode_all())


def _attributes(cloud):
    return [e for e in cloud.extensions if isinstance(e, spz.SpzExtensionCustomAttributeNiantic)]


def _assert_same_core(a, b):
    for field in ("positions", "scales", "rotations", "alphas", "colors", "sh"):
        np.testing.assert_array_equal(getattr(a, field), getattr(b, field), err_msg=field)


# ---- File layout helpers ----


def _parse(data):
    magic, version, num_points, sh_degree, frac, flags, num_streams, tbo = HEADER.unpack_from(data)
    toc = [struct.unpack_from("<QQ", data, tbo + 16 * i) for i in range(num_streams)]
    offset = tbo + 16 * num_streams
    starts = []
    for compressed, _ in toc:
        starts.append(offset)
        offset += compressed
    return {"flags": flags, "num_streams": num_streams, "tbo": tbo, "toc": toc, "starts": starts}


def _records(data, tbo):
    records, pos = [], 32
    while pos < tbo:
        rtype, length = struct.unpack_from("<II", data, pos)
        records.append((rtype, pos + 8, length))
        pos += 8 + length
    return records


def _frame_table_offset(data, payload_start):
    pos = payload_start
    (name_length,) = struct.unpack_from("<H", data, pos)
    pos += 2 + name_length
    dim, _, bits = struct.unpack_from("<HBB", data, pos)
    pos += 4
    if 1 <= bits <= 16:
        pos += 2 * 4 * dim
    return pos  # points at u32 frameCount


def _raw_zstd_frame(payload):
    """A valid zstd frame built from raw (uncompressed) blocks, so tests need no zstd library."""
    out = bytearray(struct.pack("<IBI", ZSTD_MAGIC, 0xA0, len(payload)))  # single segment, 4-byte size
    blocks = [payload[i:i + MAX_RAW_BLOCK] for i in range(0, len(payload), MAX_RAW_BLOCK)] or [b""]
    for i, block in enumerate(blocks):
        header = (len(block) << 3) | (1 if i == len(blocks) - 1 else 0)
        out += header.to_bytes(3, "little") + block
    return bytes(out)


def _rle_zstd_frame(byte, size):
    """A valid zstd frame of `size` copies of `byte`, a few bytes per 128 KiB block."""
    out = bytearray(struct.pack("<IBI", ZSTD_MAGIC, 0xA0, size))
    remaining = size
    while remaining:
        block = min(remaining, MAX_RAW_BLOCK)
        remaining -= block
        header = (block << 3) | (1 << 1) | (0 if remaining else 1)  # RLE block
        out += header.to_bytes(3, "little") + bytes([byte])
    return bytes(out)


def _inject_attribute(plain, name, dim, chunks, bits=32, ranges=None, semantic=CUSTOM,
                      frame_magic=SKIPPABLE_MAGIC, offset_delta=0, uncompressed_delta=0,
                      frame_count=None, trailing=b"", zframes=None):
    """Adds a custom attribute record and its frames to a v4 file that has no extensions.

    name may be bytes; dim and semantic are stored as given. chunks are the uncompressed payload
    pieces, each wrapped in a raw zstd frame unless zframes supplies the frames; trailing is
    appended after the frame table; the remaining knobs corrupt the result for negative tests.
    """
    info = _parse(plain)
    assert info["tbo"] == 32 and info["flags"] & FLAG_HAS_EXTENSIONS == 0
    host_compressed, _ = info["toc"][-1]

    frames, locations, offset = b"", [], host_compressed
    for i, chunk in enumerate(chunks):
        zframe = zframes[i] if zframes else _raw_zstd_frame(chunk)
        locations.append((offset + offset_delta, len(zframe), len(chunk) + uncompressed_delta))
        frames += struct.pack("<II", frame_magic, len(zframe)) + zframe
        offset += 8 + len(zframe)

    name_bytes = name if isinstance(name, bytes) else name.encode()
    payload = (struct.pack("<H", len(name_bytes)) + name_bytes
               + struct.pack("<HBB", dim, semantic, bits))
    if ranges is not None:
        n = len(ranges[0])
        payload += struct.pack(f"<{n}f", *ranges[0]) + struct.pack(f"<{n}f", *ranges[1])
    payload += struct.pack("<I", len(locations) if frame_count is None else frame_count)
    for loc in locations:
        payload += struct.pack("<QQQ", *loc)
    payload += trailing
    record = struct.pack("<II", CUSTOM_ATTRIBUTE_TYPE, len(payload)) + payload

    header = bytearray(plain[:32])
    magic, version, num_points, sh_degree, frac, flags, num_streams, _ = HEADER.unpack_from(header)
    HEADER.pack_into(header, 0, magic, version, num_points, sh_degree, frac,
                     flags | FLAG_HAS_EXTENSIONS, num_streams, 32 + len(record))
    toc = bytearray(plain[32:32 + 16 * num_streams])
    struct.pack_into("<Q", toc, 16 * (num_streams - 1), host_compressed + len(frames))
    return bytes(header) + record + bytes(toc) + plain[32 + 16 * num_streams:] + frames


def _float_chunks(values, components, num_points, pieces):
    planar = np.ascontiguousarray(
        np.asarray(values, np.float32).reshape(num_points, components).T).tobytes()
    step = -(-len(planar) // pieces)
    return [planar[i:i + step] for i in range(0, len(planar), step)]


# ---- Round trips ----


@pytest.mark.parametrize("dim", [1, 2, 3, 7, 64, 768])
def test_float_round_trip_is_exact(dim):
    cloud = _make_cloud()
    values = np.random.default_rng(1).normal(size=cloud.num_points * dim).astype(np.float32)
    cloud.extensions = [_attribute("feature", dim, values)]
    _, loaded = _round_trip(cloud)
    (attr,) = _attributes(loaded)
    assert (attr.name, attr.dim, attr.bits) == ("feature", dim, 32)
    np.testing.assert_array_equal(attr.values, values)


@pytest.mark.parametrize("bits", range(1, 17))
def test_quantized_round_trip_within_half_step(bits):
    cloud = _make_cloud()
    values = np.random.default_rng(2).uniform(-3, 5, size=cloud.num_points * 2).astype(np.float32)
    cloud.extensions = [_attribute("_Q", 2, values, bits=bits)]
    _, loaded = _round_trip(cloud)
    (attr,) = _attributes(loaded)
    per_component = values.reshape(-1, 2)
    half_step = (per_component.max(0) - per_component.min(0)) / ((1 << bits) - 1) / 2
    err = np.abs(attr.values.reshape(-1, 2) - per_component)
    assert np.all(err <= half_step * 1.0001 + 1e-6)


def test_constant_quantized_component_round_trips():
    cloud = _make_cloud()
    values = np.full(cloud.num_points, 4.25, np.float32)
    cloud.extensions = [_attribute("_CONST", 1, values, bits=4)]
    _, loaded = _round_trip(cloud)
    np.testing.assert_array_equal(_attributes(loaded)[0].values, values)


def test_multiple_attributes_and_other_extensions_coexist():
    cloud = _make_cloud()
    n = cloud.num_points
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 2.5
    a = np.arange(n, dtype=np.float32)
    b = np.linspace(0, 1, n * 3, dtype=np.float32)
    cloud.extensions = [_attribute("_A", 1, a), orbit, _attribute("_B", 3, b, bits=8)]
    _, loaded = _round_trip(cloud)
    attrs = {attr.name: attr for attr in _attributes(loaded)}
    np.testing.assert_array_equal(attrs["_A"].values, a)
    np.testing.assert_allclose(attrs["_B"].values, b, atol=1 / 255)
    orbits = [e for e in loaded.extensions if isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe)]
    assert len(orbits) == 1 and orbits[0].safe_orbit_radius_min == pytest.approx(2.5)


def test_quantized_embedding_round_trips():
    cloud = _make_cloud()
    dim = 512
    values = np.random.default_rng(4).normal(size=cloud.num_points * dim).astype(np.float32)
    cloud.extensions = [_attribute("embedding", dim, values, bits=8)]
    _, loaded = _round_trip(cloud)
    (attr,) = _attributes(loaded)
    per_component = values.reshape(-1, dim)
    half_step = (per_component.max(0) - per_component.min(0)) / 255 / 2
    assert np.all(np.abs(attr.values.reshape(-1, dim) - per_component) <= half_step * 1.0001 + 1e-6)


@pytest.mark.parametrize("name", ["roughness", "nx", "f_rest_0", "_TEMPERATURE", "my attr #1"])
def test_names_are_free_form_labels(name):
    cloud = _make_cloud()
    values = np.arange(cloud.num_points, dtype=np.float32)
    cloud.extensions = [_attribute(name, 1, values)]
    _, loaded = _round_trip(cloud)
    (attr,) = _attributes(loaded)
    assert (attr.name, attr.semantic) == (name, S.CUSTOM)
    np.testing.assert_array_equal(attr.values, values)


def test_sh_degree_zero_hosts_frames_in_rotations_stream():
    cloud = _make_cloud(sh_degree=0)
    values = np.arange(cloud.num_points, dtype=np.float32)
    cloud.extensions = [_attribute("_ID", 1, values)]
    data, loaded = _round_trip(cloud)
    assert _parse(data)["num_streams"] == 5
    np.testing.assert_array_equal(_attributes(loaded)[0].values, values)


def test_core_data_is_unchanged_by_attributes():
    cloud = _make_cloud()
    _, plain = _round_trip(cloud)
    cloud.extensions = [_attribute("_X", 1, np.ones(cloud.num_points, np.float32))]
    _, with_attr = _round_trip(cloud)
    _assert_same_core(plain, with_attr)


# ---- File layout ----


def test_frames_are_appended_to_last_core_stream():
    cloud = _make_cloud()
    plain = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    cloud.extensions = [_attribute("_X", 1, np.arange(cloud.num_points, dtype=np.float32))]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())

    p, info = _parse(plain), _parse(data)
    assert info["num_streams"] == p["num_streams"] == 6
    assert info["flags"] & FLAG_HAS_EXTENSIONS
    # Every core TOC entry is unchanged except the host's compressed size, which grows by the frames.
    assert info["toc"][:-1] == p["toc"][:-1]
    assert info["toc"][-1][1] == p["toc"][-1][1]
    assert info["toc"][-1][0] > p["toc"][-1][0]

    ((rtype, payload_start, _),) = _records(data, info["tbo"])
    assert rtype == CUSTOM_ATTRIBUTE_TYPE
    table = _frame_table_offset(data, payload_start)
    (count,) = struct.unpack_from("<I", data, table)
    assert count == 1
    frame_offset, compressed, uncompressed = struct.unpack_from("<QQQ", data, table + 4)
    assert frame_offset == p["toc"][-1][0]  # right after the core frame
    assert uncompressed == cloud.num_points * 4
    frame = info["starts"][-1] + frame_offset
    assert struct.unpack_from("<II", data, frame) == (SKIPPABLE_MAGIC, compressed)
    assert struct.unpack_from("<I", data, frame + 8)[0] == ZSTD_MAGIC


# ---- Coordinate conversion ----


def test_spatial_converts_like_positions_and_custom_never_converts():
    cloud = _make_cloud()
    directions = cloud.positions.copy()
    cloud.extensions = [_attribute("normal", 3, directions, semantic=S.SPATIAL),
                        _attribute("velocity", 3, directions, semantic=S.SPATIAL),
                        _attribute("raw", 3, directions)]
    pack = spz.PackOptions()
    pack.from_coord = spz.RDF
    unpack = _decode_all()
    unpack.to_coord = spz.RUB
    _, loaded = _round_trip(cloud, pack, unpack)

    attrs = {attr.name: attr for attr in _attributes(loaded)}
    # Positions are stored as 12-bit fixed point; the attributes are lossless.
    for name in ("normal", "velocity"):
        np.testing.assert_allclose(attrs[name].values, loaded.positions, atol=1 / 4096)
        assert not np.allclose(attrs[name].values, directions)
    np.testing.assert_array_equal(attrs["raw"].values, directions)
    assert attrs["normal"].semantic == S.SPATIAL and attrs["raw"].semantic == S.CUSTOM
    # The caller's cloud is left alone.
    np.testing.assert_array_equal(_attributes(cloud)[0].values, directions)


def test_spatial_converts_every_triple():
    cloud = _make_cloud()
    two_vectors = np.concatenate([cloud.positions.reshape(-1, 3)] * 2, axis=1).reshape(-1)
    cloud.extensions = [_attribute("frame", 6, two_vectors, semantic=S.SPATIAL)]
    cloud.convert_coordinates(spz.RUB, spz.RDF)
    expected = np.concatenate([cloud.positions.reshape(-1, 3)] * 2, axis=1).reshape(-1)
    np.testing.assert_array_equal(_attributes(cloud)[0].values, expected)


def test_spatial_follows_convert_coordinates():
    cloud = _make_cloud()
    directions = cloud.positions.copy()
    cloud.extensions = [_attribute("normal", 3, directions, semantic=S.SPATIAL)]
    cloud.convert_coordinates(spz.RUB, spz.RDF)
    np.testing.assert_array_equal(_attributes(cloud)[0].values, cloud.positions)


# ---- Save-time validation ----


@pytest.mark.parametrize(
    "name,semantic,dim,bits,values_per_point",
    [
        ("", S.CUSTOM, 1, 32, 1),  # empty name
        ("Gr\u00f6\u00dfe", S.CUSTOM, 1, 32, 1),  # not printable ASCII
        ("tab\there", S.CUSTOM, 1, 32, 1),  # control character
        ("x", S.CUSTOM, 0, 32, 0),  # dim 0
        ("n", S.SPATIAL, 4, 32, 4),  # SPATIAL needs a multiple of 3
        ("v", S.SPATIAL, 1, 32, 1),  # SPATIAL needs a multiple of 3
        ("x", S.CUSTOM, 1, 17, 1),  # bad bits
        ("x", S.CUSTOM, 1, 31, 1),  # bad bits
        ("x", S.CUSTOM, 1, 0, 1),  # bad bits
        ("x", S.CUSTOM, 2, 32, 1),  # wrong value count
    ],
)
def test_invalid_attribute_fails_save(name, semantic, dim, bits, values_per_point):
    cloud = _make_cloud()
    values = np.zeros(cloud.num_points * values_per_point)
    cloud.extensions = [_attribute(name, dim, values, bits, semantic)]
    with pytest.raises(RuntimeError):
        spz.save_spz_to_buffer(cloud, spz.PackOptions())


def test_duplicate_names_fail_save():
    cloud = _make_cloud()
    n = cloud.num_points
    cloud.extensions = [_attribute("_X", 1, np.zeros(n)), _attribute("_X", 1, np.ones(n))]
    with pytest.raises(RuntimeError):
        spz.save_spz_to_buffer(cloud, spz.PackOptions())


def test_non_finite_quantized_values_fail_save():
    cloud = _make_cloud()
    values = np.zeros(cloud.num_points, np.float32)
    values[3] = np.nan
    cloud.extensions = [_attribute("_X", 1, values, bits=8)]
    with pytest.raises(RuntimeError):
        spz.save_spz_to_buffer(cloud, spz.PackOptions())


def test_legacy_version_fails_save():
    cloud = _make_cloud()
    cloud.extensions = [_attribute("_X", 1, np.zeros(cloud.num_points))]
    pack = spz.PackOptions()
    pack.version = 3
    with pytest.raises(RuntimeError):
        spz.save_spz_to_buffer(cloud, pack)


# ---- Reading hand-built and malformed files ----


@pytest.fixture
def plain_file():
    cloud = _make_cloud()
    return cloud, spz.save_spz_to_buffer(cloud, spz.PackOptions())


@pytest.mark.parametrize("pieces", [1, 2, 3])
def test_reads_attribute_split_across_frames(plain_file, pieces):
    cloud, plain = plain_file
    values = np.random.default_rng(3).normal(size=cloud.num_points * 3).astype(np.float32)
    data = _inject_attribute(plain, "_SPLIT", 3, _float_chunks(values, 3, cloud.num_points, pieces))
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    np.testing.assert_array_equal(_attributes(loaded)[0].values, values)


def test_reads_quantized_attribute_with_stored_ranges(plain_file):
    cloud, plain = plain_file
    levels = (np.arange(cloud.num_points) % 16).astype(np.uint8)
    data = _inject_attribute(plain, "_Q", 1, [levels.tobytes()], bits=4, ranges=([-1.0], [2.0]))
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    np.testing.assert_allclose(_attributes(loaded)[0].values, -1.0 + levels * (3.0 / 15), atol=1e-6)


@pytest.mark.parametrize(
    "corruption",
    [
        {"frame_magic": 0x184D2A50},  # a different skippable magic
        {"offset_delta": 1},  # offset doesn't point at the frame header
        {"offset_delta": -8},  # offset points into the core frame
        {"offset_delta": 1 << 40},  # offset past the host stream
        {"uncompressed_delta": 4},  # claims more data than the point count needs
        {"frame_count": 5},  # frame table truncated
    ],
)
def test_malformed_attribute_is_dropped_and_core_still_loads(plain_file, corruption):
    cloud, plain = plain_file
    values = np.arange(cloud.num_points, dtype=np.float32)
    data = _inject_attribute(plain, "_BAD", 1, _float_chunks(values, 1, cloud.num_points, 1),
                             **corruption)
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    _assert_same_core(loaded, spz.load_spz_from_buffer(plain, _decode_all()))


@pytest.mark.parametrize(
    "semantic,dim",
    [
        (2, 1),  # a semantic this reader doesn't know: its conversion rules are unknown
        (SPATIAL, 1),  # SPATIAL whose dim isn't a multiple of 3
        (CUSTOM, 0),  # dim 0
    ],
)
def test_attribute_with_invalid_descriptor_is_dropped(plain_file, semantic, dim):
    cloud, plain = plain_file
    chunk = np.zeros(cloud.num_points * max(dim, 1), np.float32).tobytes()
    data = _inject_attribute(plain, "x", dim, [chunk], semantic=semantic)
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def _byte_planes(levels):
    """One component's 16-bit levels as stored: every low byte, then every high byte."""
    levels = np.asarray(levels, dtype=np.uint16)
    return (levels & 0xFF).astype(np.uint8).tobytes() + (levels >> 8).astype(np.uint8).tobytes()


def test_reads_16_bit_levels_from_byte_planes(plain_file):
    cloud, plain = plain_file
    levels = (np.arange(cloud.num_points) * 1297) % 65536  # exercises both bytes of every level
    data = _inject_attribute(plain, "_W", 1, [_byte_planes(levels)], bits=16,
                             ranges=([-2.0], [2.0]))
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    np.testing.assert_allclose(_attributes(loaded)[0].values, -2.0 + levels * (4.0 / 65535),
                               atol=1e-6)


def test_writer_matches_the_byte_plane_reader():
    # The reader is pinned to the spec above, so exact levels surviving a round trip pin the writer.
    cloud = _make_cloud()
    levels = (np.arange(cloud.num_points) * 1297) % 65536
    levels[0], levels[1] = 0, 65535  # the stored range is then exactly [-2, 2]
    expected = -2.0 + levels * (4.0 / 65535)
    cloud.extensions = [_attribute("_W", 1, expected.astype(np.float32), bits=16)]
    data, loaded = _round_trip(cloud)
    np.testing.assert_allclose(_attributes(loaded)[0].values, expected, atol=1e-6)
    ((_, payload_start, _),) = _records(data, _parse(data)["tbo"])
    table = _frame_table_offset(data, payload_start)
    (frame_count,) = struct.unpack_from("<I", data, table)
    sizes = [struct.unpack_from("<QQQ", data, table + 4 + 24 * i)[2] for i in range(frame_count)]
    assert sum(sizes) == cloud.num_points * 2


def test_twelve_bit_level_above_range_drops_attribute(plain_file):
    cloud, plain = plain_file
    levels = np.full(cloud.num_points, 5000)  # 12-bit attribute allows at most 4095
    data = _inject_attribute(plain, "_Q", 1, [_byte_planes(levels)], bits=12,
                             ranges=([0.0], [1.0]))
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_16_bit_attribute_with_one_byte_per_level_is_dropped(plain_file):
    cloud, plain = plain_file
    levels = np.zeros(cloud.num_points, np.uint8)  # half the bytes a 16-bit attribute needs
    data = _inject_attribute(plain, "_Q", 1, [levels.tobytes()], bits=16, ranges=([0.0], [1.0]))
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_out_of_range_quantized_level_drops_attribute(plain_file):
    cloud, plain = plain_file
    levels = np.full(cloud.num_points, 200, np.uint8)  # 4-bit attribute allows at most 15
    data = _inject_attribute(plain, "_Q", 1, [levels.tobytes()], bits=4, ranges=([0.0], [1.0]))
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def _two_attribute_file(cloud):
    cloud.extensions = [_attribute("_A", 1, np.ones(cloud.num_points)),
                        _attribute("_B", 1, np.zeros(cloud.num_points))]
    data = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    (_, start_a, _), (_, start_b, _) = _records(data, _parse(data)["tbo"])
    return bytearray(data), start_a, start_b


def test_duplicate_names_in_file_drop_both_attributes(plain_file):
    cloud, _ = plain_file
    data, _, start_b = _two_attribute_file(cloud)
    data[start_b + 2:start_b + 4] = b"_A"  # rename _B to _A (same length)
    loaded = spz.load_spz_from_buffer(bytes(data), _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_non_ascii_name_in_file_drops_attribute(plain_file):
    cloud, plain = plain_file
    chunk = np.zeros(cloud.num_points, np.float32).tobytes()
    data = _inject_attribute(plain, b"\xff\xfe", 1, [chunk])
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_bytes_after_frame_table_are_ignored(plain_file):
    cloud, plain = plain_file
    values = np.arange(cloud.num_points, dtype=np.float32)
    data = _inject_attribute(plain, "x", 1, [values.tobytes()], trailing=b"future fields")
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    np.testing.assert_array_equal(_attributes(loaded)[0].values, values)


def test_highly_compressed_frame_within_size_cap_loads(plain_file):
    cloud, plain = plain_file
    dim = 1000  # 200 KB of zeros in a few hundred bytes
    size = cloud.num_points * dim * 4
    data = _inject_attribute(plain, "zeros", dim, [b"\0" * size], zframes=[_rle_zstd_frame(0, size)])
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    np.testing.assert_array_equal(_attributes(loaded)[0].values, np.zeros(cloud.num_points * dim))


def test_decoded_size_far_beyond_file_size_is_dropped_before_allocating(plain_file):
    cloud, plain = plain_file
    dim = 65535  # 13 MB claimed by a frame of a few KB: over 1024x the file size
    size = cloud.num_points * dim * 4
    data = _inject_attribute(plain, "huge", dim, [b""], zframes=[_rle_zstd_frame(0, size)],
                             uncompressed_delta=size)
    assert size > 1024 * len(data)
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_frame_reused_by_many_references_is_capped_by_file_size(plain_file):
    cloud, plain = plain_file
    host_compressed, _ = _parse(plain)["toc"][-1]
    per_frame = cloud.num_points * 4 * 300
    references = 200
    zframe = _rle_zstd_frame(0, per_frame)
    location = struct.pack("<QQQ", host_compressed, len(zframe), per_frame)
    data = _inject_attribute(plain, "reused", 300 * references, [b""], zframes=[zframe],
                             uncompressed_delta=per_frame, frame_count=references,
                             trailing=location * (references - 1))
    assert per_frame * references > 1024 * len(data)
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_find_custom_attribute():
    cloud = _make_cloud()
    values = np.arange(cloud.num_points, dtype=np.float32)
    cloud.extensions = [_attribute("a", 1, values), _attribute("b", 1, values * 2)]
    _, loaded = _round_trip(cloud)
    np.testing.assert_array_equal(spz.find_custom_attribute(loaded, "b").values, values * 2)
    assert spz.find_custom_attribute(loaded, "missing") is None


def test_values_accept_two_dimensional_arrays():
    cloud = _make_cloud()
    attr = spz.SpzExtensionCustomAttributeNiantic()
    attr.name, attr.dim = "normal", 3
    attr.values = cloud.positions.reshape(-1, 3)
    np.testing.assert_array_equal(attr.values, cloud.positions)


def test_duplicate_name_drops_both_even_when_one_is_invalid(plain_file):
    cloud, _ = plain_file
    data, start_a, start_b = _two_attribute_file(cloud)
    data[start_b + 2:start_b + 4] = b"_A"  # rename _B to _A
    table_b = _frame_table_offset(data, start_b)
    offset = struct.unpack_from("<Q", data, table_b + 4)[0]
    struct.pack_into("<Q", data, table_b + 4, offset + 1)  # and break its frame
    loaded = spz.load_spz_from_buffer(bytes(data), _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_payload_must_be_exactly_one_zstd_frame(plain_file):
    cloud, plain = plain_file
    values = np.arange(cloud.num_points, dtype=np.float32).tobytes()
    padded = _raw_zstd_frame(values) + struct.pack("<II", SKIPPABLE_MAGIC, 0)
    data = _inject_attribute(plain, "x", 1, [values], zframes=[padded])
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    assert _attributes(loaded) == []
    assert loaded.num_points == cloud.num_points


def test_duplicate_with_malformed_descriptor_does_not_count(plain_file):
    cloud, _ = plain_file
    data, _, start_b = _two_attribute_file(cloud)
    data[start_b + 2:start_b + 4] = b"_A"  # rename _B to _A
    struct.pack_into("<H", data, start_b + 4, 0)  # and give it dim 0, a malformed descriptor
    loaded = spz.load_spz_from_buffer(bytes(data), _decode_all())
    assert [attr.name for attr in _attributes(loaded)] == ["_A"]


def test_attribute_larger_than_one_frame_is_split_and_round_trips():
    cloud = _make_cloud(num_points=5000)
    dim = 1000  # 20 MB of float32: more than one 16 MiB frame
    values = np.random.default_rng(5).normal(size=cloud.num_points * dim).astype(np.float32)
    cloud.extensions = [_attribute("big", dim, values)]
    data, loaded = _round_trip(cloud)
    ((_, payload_start, _),) = _records(data, _parse(data)["tbo"])
    table = _frame_table_offset(data, payload_start)
    (frame_count,) = struct.unpack_from("<I", data, table)
    sizes = [struct.unpack_from("<QQQ", data, table + 4 + 24 * i)[2] for i in range(frame_count)]
    assert frame_count == 2 and sizes == [16 << 20, values.nbytes - (16 << 20)]
    np.testing.assert_array_equal(_attributes(loaded)[0].values, values)


def test_reads_extension_frame_larger_than_this_writer_produces():
    # This writer caps extension frames at 16 MiB; readers must accept any size within the format.
    cloud = _make_cloud(num_points=5000)
    plain = spz.save_spz_to_buffer(cloud, spz.PackOptions())
    dim = 1000
    size = cloud.num_points * dim * 4  # 20 MB in one extension frame, a few KB on disk
    data = _inject_attribute(plain, "big", dim, [b""], zframes=[_rle_zstd_frame(0, size)],
                             uncompressed_delta=size)
    assert size > 16 << 20
    loaded = spz.load_spz_from_buffer(data, _decode_all())
    (attr,) = _attributes(loaded)
    np.testing.assert_array_equal(attr.values, np.zeros(cloud.num_points * dim, np.float32))


# ---- Choosing which attributes to decode ----


def _unpack_only(names):
    options = spz.UnpackOptions()
    options.custom_attributes = names
    return options


def _cloud_with_attributes_and_orbit():
    cloud = _make_cloud()
    orbit = spz.SpzExtensionSafeOrbitCameraAdobe()
    orbit.safe_orbit_radius_min = 2.5
    a = np.arange(cloud.num_points, dtype=np.float32)
    b = np.linspace(0, 1, cloud.num_points * 3, dtype=np.float32)
    cloud.extensions = [_attribute("_A", 1, a), orbit, _attribute("_B", 3, b)]
    return cloud, a, b


def test_none_decodes_every_attribute():
    cloud, _, _ = _cloud_with_attributes_and_orbit()
    _, loaded = _round_trip(cloud, unpack=_unpack_only(None))
    assert sorted(attr.name for attr in _attributes(loaded)) == ["_A", "_B"]


@pytest.mark.parametrize("unpack", [spz.UnpackOptions, lambda: _unpack_only([])])
def test_decodes_no_attributes_by_default_and_keeps_everything_else(unpack):
    assert spz.UnpackOptions().custom_attributes == []
    cloud, _, _ = _cloud_with_attributes_and_orbit()
    data, loaded = _round_trip(cloud, unpack=unpack())
    assert _attributes(loaded) == []
    orbits = [e for e in loaded.extensions if isinstance(e, spz.SpzExtensionSafeOrbitCameraAdobe)]
    assert len(orbits) == 1 and orbits[0].safe_orbit_radius_min == pytest.approx(2.5)
    _assert_same_core(loaded, spz.load_spz_from_buffer(data, _decode_all()))


def test_decodes_only_the_listed_attributes():
    cloud, _, b = _cloud_with_attributes_and_orbit()
    _, loaded = _round_trip(cloud, unpack=_unpack_only(["_B", "_MISSING"]))
    (attr,) = _attributes(loaded)
    assert attr.name == "_B"
    np.testing.assert_array_equal(attr.values, b)


def test_selection_applies_to_load_spz_from_file(tmp_path):
    cloud, a, _ = _cloud_with_attributes_and_orbit()
    path = str(tmp_path / "attributes.spz")
    assert spz.save_spz(cloud, spz.PackOptions(), path)
    (attr,) = _attributes(spz.load_spz(path, _unpack_only(["_A"])))
    np.testing.assert_array_equal(attr.values, a)


def test_unlisted_attribute_frames_are_never_checked_or_decompressed(plain_file, capfd):
    cloud, plain = plain_file
    values = np.arange(cloud.num_points, dtype=np.float32)
    data = _inject_attribute(plain, "_BAD", 1, _float_chunks(values, 1, cloud.num_points, 1),
                             offset_delta=1 << 40)
    capfd.readouterr()
    spz.load_spz_from_buffer(data, _decode_all())
    assert "Dropped extension" in capfd.readouterr().out
    loaded = spz.load_spz_from_buffer(data, _unpack_only([]))
    assert "Dropped extension" not in capfd.readouterr().out
    _assert_same_core(loaded, spz.load_spz_from_buffer(plain, _decode_all()))


def test_requested_duplicate_names_are_still_dropped(plain_file):
    cloud, _ = plain_file
    data, _, start_b = _two_attribute_file(cloud)
    data[start_b + 2:start_b + 4] = b"_A"  # rename _B to _A (same length)
    loaded = spz.load_spz_from_buffer(bytes(data), _unpack_only(["_A"]))
    assert _attributes(loaded) == []
