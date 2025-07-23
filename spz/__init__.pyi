"""Type stubs for SPZ Python bindings."""

from typing import List, overload

# Coordinate System Enum
class CoordinateSystem:
    UNSPECIFIED: int
    LDB: int  # Left Down Back
    RDB: int  # Right Down Back
    LUB: int  # Left Up Back
    RUB: int  # Right Up Back, Three.js coordinate system
    LDF: int  # Left Down Front
    RDF: int  # Right Down Front, PLY coordinate system
    LUF: int  # Left Up Front, GLB coordinate system
    RUF: int  # Right Up Front, Unity coordinate system

# Module-level constants
UNSPECIFIED: int
LDB: int
RDB: int
LUB: int
RUB: int
LDF: int
RDF: int
LUF: int
RUF: int

class CoordinateConverter:
    def __init__(self) -> None: ...
    flipP: bool
    flipQ: bool
    flipSh: bool

class PackOptions:
    def __init__(self) -> None: ...
    from_: int  # Note: 'from' is a Python keyword in the C++ binding
    sh1Bits: int  # Bits for SH degree 1 coefficients
    shRestBits: int  # Bits for SH degree 2+ coefficients
    hasSafeOrbit: bool  # Whether safe orbit data is present
    safeOrbitElevationMin: float  # Minimum elevation for safe orbit (radians)
    safeOrbitElevationMax: float  # Maximum elevation for safe orbit (radians)
    safeOrbitRadiusMin: float  # Minimum radius for safe orbit

class UnpackOptions:
    def __init__(self) -> None: ...
    to: int

class UnpackedGaussian:
    def __init__(self) -> None: ...
    position: List[float]  # 3D vector
    rotation: List[float]  # quaternion/4D vector
    scale: List[float]  # 3D vector
    color: List[float]  # RGB color
    alpha: float
    shR: List[float]  # Spherical harmonics for Red
    shG: List[float]  # Spherical harmonics for Green
    shB: List[float]  # Spherical harmonics for Blue

class PackedGaussian:
    def __init__(self) -> None: ...
    def unpack(self, usesFloat16: bool, fractionalBits: int, c: CoordinateConverter, shMin: float = -1.0, shMax: float = 1.0) -> UnpackedGaussian: ...

class PackedGaussians:
    def __init__(self) -> None: ...
    numPoints: int
    shDegree: int
    fractionalBits: int
    antialiased: bool
    sh1Bits: int  # Bits for SH degree 1 coefficients
    shRestBits: int  # Bits for SH degree 2+ coefficients
    shMin: float  # Minimum SH coefficient value used for quantization
    shMax: float  # Maximum SH coefficient value used for quantization
    hasSafeOrbit: bool  # Whether safe orbit data is present
    safeOrbitElevationMin: float  # Minimum elevation for safe orbit (radians)
    safeOrbitElevationMax: float  # Maximum elevation for safe orbit (radians)
    safeOrbitRadiusMin: float  # Minimum radius for safe orbit
    positions: List[float]
    scales: List[float]
    rotations: List[float]
    alphas: List[float]
    colors: List[float]
    sh: List[float]

    def usesFloat16(self) -> bool: ...
    def at(self, index: int) -> PackedGaussian: ...
    def unpack(self, i: int, c: CoordinateConverter) -> UnpackedGaussian: ...

class SpzFloatBuffer:
    def __init__(self) -> None: ...
    count: int
    data: List[float]

class GaussianCloudData:
    def __init__(self) -> None: ...
    numPoints: int
    shDegree: int
    antialiased: bool
    positions: List[float]
    scales: List[float]
    rotations: List[float]
    alphas: List[float]
    colors: List[float]
    sh: List[float]

class GaussianCloud:
    def __init__(self) -> None: ...
    numPoints: int
    shDegree: int
    antialiased: bool
    positions: List[float]
    scales: List[float]
    rotations: List[float]
    alphas: List[float]
    colors: List[float]
    sh: List[float]

    def data(self) -> GaussianCloudData: ...
    def convertCoordinates(self, from_system: int, to_system: int) -> None: ...
    def rotate180DegAboutX(self) -> None: ...
    def medianVolume(self) -> float: ...

def coordinateConverter(from_system: int, to_system: int) -> CoordinateConverter: ...

# saveSpz overloads
@overload
def saveSpz(g: GaussianCloud, options: PackOptions, output: List[int]) -> bool: ...
@overload
def saveSpz(g: GaussianCloud, options: PackOptions, filename: str) -> bool: ...
@overload
def saveSpz(g: GaussianCloud, filename: str) -> bool: ...  # Default options
@overload
def saveSpz(g: GaussianCloud, output: List[int]) -> bool: ...  # Default options

# loadSpz overloads
@overload
def loadSpz(data: List[int], options: UnpackOptions) -> GaussianCloud: ...
@overload
def loadSpz(filename: str, options: UnpackOptions) -> GaussianCloud: ...
@overload
def loadSpz(filename: str) -> GaussianCloud: ...  # Default options
@overload
def loadSpz(data: List[int]) -> GaussianCloud: ...  # Default options

# loadSpzPacked overloads
@overload
def loadSpzPacked(filename: str) -> PackedGaussians: ...
@overload
def loadSpzPacked(data: bytes, size: int) -> PackedGaussians: ...
@overload
def loadSpzPacked(data: List[int]) -> PackedGaussians: ...

# PLY functions
@overload
def saveSplatToPly(gaussians: GaussianCloud, options: PackOptions, filename: str) -> None: ...
@overload
def saveSplatToPly(gaussians: GaussianCloud, filename: str) -> None: ...  # Default options
@overload
def loadSplatFromPly(filename: str, options: UnpackOptions) -> GaussianCloud: ...
@overload
def loadSplatFromPly(filename: str) -> GaussianCloud: ...  # Default options

# Utility functions
@overload
def saveSpzToBytes(gaussians: GaussianCloud, options: PackOptions) -> bytes: ...
@overload
def saveSpzToBytes(gaussians: GaussianCloud) -> bytes: ...  # Default options
def serializePackedGaussians(packed: PackedGaussians) -> bytes: ...
def compressGzipped(data: bytes) -> bytes: ...
