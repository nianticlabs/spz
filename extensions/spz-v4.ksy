# Kaitai Struct grammar for SPZ version 4 (NGSP) files and their extension records.
#
# It describes the byte layout only: the core streams and extension frames are zstd data, and
# validity rules (e.g. when a temporal record must be ignored) are in README.md and
# extensions/README.md. Open it in the Kaitai Web IDE (https://ide.kaitai.io) with a file from
# tests/conformance to browse a file's structure. tests/python/test_grammar.py checks that its
# record layouts agree with tools/4dgs/spz_layout.py, which reads the conformance corpus.
meta:
  id: spz_v4
  title: SPZ Gaussian splat file, version 4
  file-extension: spz
  license: MIT
  endian: le

seq:
  - id: header
    type: header
  - id: extension_zone
    size: header.toc_byte_offset - 32
    type: extension_records
    if: (header.flags & 2) != 0
  - id: toc
    type: toc_entry
    repeat: expr
    repeat-expr: header.num_streams
  - id: streams
    size: toc[_index].compressed_size
    repeat: expr
    repeat-expr: header.num_streams
    doc: |
      zstd frames, in order: positions, alphas, colors, scales, rotations, SH (absent at degree 0).
      The last stream's bytes also hold the extension frames that framed records point at.

types:
  header:
    seq:
      - id: magic
        contents: "NGSP"
      - id: version
        type: u4
        valid: 4
      - id: num_points
        type: u4
      - id: sh_degree
        type: u1
      - id: fractional_bits
        type: u1
      - id: flags
        type: u1
        doc: "Bit 0: antialiased. Bit 1: extension records follow the header."
      - id: num_streams
        type: u1
      - id: toc_byte_offset
        type: u4
      - id: reserved
        size: 12

  toc_entry:
    seq:
      - id: compressed_size
        type: u8
      - id: uncompressed_size
        type: u8

  extension_records:
    seq:
      - id: records
        type: extension_record
        repeat: eos

  extension_record:
    seq:
      - id: type
        type: u4
        enum: extension_type
      - id: byte_length
        type: u4
      - id: payload
        size: byte_length
        type:
          switch-on: type
          cases:
            "extension_type::adobe_safe_orbit_camera": safe_orbit_camera
            "extension_type::adobe_coordinate_system": coordinate_system
            "extension_type::niantic_custom_attribute": custom_attribute
            "extension_type::niantic_temporal": temporal

  safe_orbit_camera:
    seq:
      - id: elevation_min
        type: f4
      - id: elevation_max
        type: f4
      - id: radius_min
        type: f4

  coordinate_system:
    seq:
      - id: coordinate_system
        type: u4

  frame_table:
    doc: The extension frames holding a framed record's data, in order.
    seq:
      - id: frame_count
        type: u4
      - id: frames
        type: frame_location
        repeat: expr
        repeat-expr: frame_count

  frame_location:
    seq:
      - id: frame_offset
        type: u8
        doc: From the start of the last core stream's data to the frame's skippable-frame magic.
      - id: compressed_size
        type: u8
      - id: uncompressed_size
        type: u8

  custom_attribute:
    seq:
      - id: name_length
        type: u2
      - id: name
        type: str
        size: name_length
        encoding: ASCII
      - id: dim
        type: u2
      - id: semantic
        type: u1
        doc: "0 CUSTOM, 1 SPATIAL."
      - id: bits
        type: u1
        doc: "32 = float32; 1-16 = quantized."
      - id: min
        type: f4
        repeat: expr
        repeat-expr: dim
        if: bits != 32
      - id: max
        type: f4
        repeat: expr
        repeat-expr: dim
        if: bits != 32
      - id: frame_table
        type: frame_table
      - id: point_range
        type: point_range
        if: not _io.eof
        doc: Absent means the attribute covers every point.

  point_range:
    seq:
      - id: point_begin
        type: u4
      - id: point_end
        type: u4

  temporal:
    doc: SPZ_NIANTIC_temporal, an animated clip. See extensions/README.md for the semantics.
    seq:
      - id: clip_start
        type: f8
      - id: clip_end
        type: f8
      - id: poster_time
        type: f8
      - id: loop
        type: u1
        doc: "0 ONCE, 1 LOOP, 2 PING_PONG."
      - id: reserved
        type: u1
      - id: nominal_fps
        type: f4
      - id: opacity_cutoff
        type: f4
      - id: channel_count
        type: u1
      - id: channels
        type: temporal_channel
        repeat: expr
        repeat-expr: channel_count
      - id: driver_count
        type: u4
      - id: keyframe_count
        type: u4
      - id: keyframe_times
        type: f8
        repeat: expr
        repeat-expr: keyframe_count
      - id: name_count
        type: u4
      - id: driver_names
        type: name
        repeat: expr
        repeat-expr: name_count
      - id: frame_table
        type: frame_table
        doc: |
          Data: driver_count * 7 * keyframe_count float32 driver values, then each channel's
          data_bytes, in channel order.
      - id: time_index
        type: time_index
        if: not _io.eof
        doc: Absent in files written before it existed.

  time_index:
    doc: |
      The clip divided into bin_count time bins [boundaries[b], boundaries[b + 1]); each lists the
      points that may be visible (at least opacity_cutoff opaque) during it, a superset.
    seq:
      - id: bin_count
        type: u4
      - id: boundaries
        type: f8
        repeat: expr
        repeat-expr: bin_count + 1
      - id: bins
        type: point_range
        repeat: expr
        repeat-expr: bin_count

  temporal_channel:
    seq:
      - id: descriptor_bytes
        type: u4
      - id: descriptor
        size: descriptor_bytes
        type: channel_descriptor

  channel_descriptor:
    seq:
      - id: target
        type: u1
        enum: temporal_target
      - id: type
        type: u1
        enum: channel_type
      - id: flags
        type: u1
        doc: "Bit 0: REQUIRED."
      - id: dim
        type: u2
      - id: point_begin
        type: u4
      - id: point_end
        type: u4
        doc: 0 means the channel covers every point.
      - id: data_bytes
        type: u8
      - id: encoding
        type: u1
        enum: channel_encoding
      - id: bits
        type: u1
        if: encoding == channel_encoding::range
      - id: min
        type: f4
        repeat: expr
        repeat-expr: dim
        if: encoding == channel_encoding::range
      - id: max
        type: f4
        repeat: expr
        repeat-expr: dim
        if: encoding == channel_encoding::range
      - id: steps
        type: f4
        repeat: expr
        repeat-expr: dim
        if: encoding == channel_encoding::step

  name:
    seq:
      - id: length
        type: u2
      - id: value
        type: str
        size: length
        encoding: ASCII

enums:
  extension_type:
    0xadbe0002: adobe_safe_orbit_camera
    0xadbe0003: adobe_coordinate_system
    0x4e530001: niantic_georeference
    0x4e530002: niantic_custom_attribute
    0x4e530003: niantic_temporal

  temporal_target:
    1: position
    2: rotation
    3: opacity
    6: transform

  channel_type:
    1: polynomial
    2: window
    4: skin

  channel_encoding:
    0: float32
    1: range
    2: step
