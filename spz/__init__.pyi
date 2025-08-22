"""Type stubs for SPZ Python bindings."""

from abc import ABC, abstractmethod
from typing import IO, Iterable, List, Optional, overload

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
LATEST_SPZ_HEADER_VERSION: int

class SpzExtensionType:
    SPZ_ADOBE_sh_quantization: int
    SPZ_ADOBE_safe_orbit_camera: int

SPZ_ADOBE_sh_quantization: int
SPZ_ADOBE_safe_orbit_camera: int

class SpzExtensionBase(ABC):
    extensionType: SpzExtensionType

    def __init__(self, t: SpzExtensionType) -> None: ...
    @abstractmethod
    def write(self, os: IO[bytes]) -> None: ...
    @abstractmethod
    def copyAsRawData(self) -> SpzExtensionBase: ...

    @classmethod
    @abstractmethod
    def read(cls, is_: IO[bytes]) -> Optional[SpzExtensionBase]: ...
    @classmethod
    @abstractmethod
    def type(cls) -> SpzExtensionType: ...

class SpzExtensionSHQuantizationAdobe(SpzExtensionBase):
    sh1Bits: int
    shRestBits: int
    shMin: float
    shMax: float

    def __init__(self) -> None: ...
    def write(self, os: IO[bytes]) -> None: ...
    def copyAsRawData(self) -> SpzExtensionBase: ...

    @classmethod
    def read(cls, is_: IO[bytes]) -> Optional[SpzExtensionBase]: ...
    @classmethod
    def type(cls) -> SpzExtensionType: ...

class SpzExtensionSafeOrbitCameraAdobe(SpzExtensionBase):
    safeOrbitElevationMin: float
    safeOrbitElevationMax: float
    safeOrbitRadiusMin: float

    def __init__(self) -> None: ...
    def write(self, os: IO[bytes]) -> None: ...
    def copyAsRawData(self) -> SpzExtensionBase: ...

    @classmethod
    def read(cls, is_: IO[bytes]) -> Optional[SpzExtensionBase]: ...
    @classmethod
    def type(cls) -> SpzExtensionType: ...

class CoordinateConverter:
    def __init__(self) -> None: ...
    flipP: bool
    flipQ: bool
    flipSh: bool

class PackOptions:
    def __init__(self) -> None: ...
    version: int  # SPZ file format version
    from_: int  # Note: 'from' is a Python keyword in the C++ binding
    sh1Bits: int  # Bits for SH degree 1 coefficients
    shRestBits: int  # Bits for SH degree 2+ coefficients
    enableSHMinMaxScaling: bool  # Enable SH min/max scaling for quantization

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
    extensions: Iterable[SpzExtensionBase]
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

class SpzExtensionNode:
    type: SpzExtensionType
    data: object
    next: Optional[SpzExtensionNode]

    def __init__(self) -> None: ...

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
    extensions: Optional[SpzExtensionNode]

class GaussianCloud:
    def __init__(self) -> None: ...
    numPoints: int
    shDegree: int
    antialiased: bool
    extensions: Iterable[SpzExtensionBase]
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
def readAllExtensions(is_: IO[bytes]) -> List[SpzExtensionBase]: ...
def writeAllExtensions(list_: Iterable[SpzExtensionBase], os: IO[bytes]) -> None: ...

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
