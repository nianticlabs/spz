

export enum CoordinateSystem {
  UNSPECIFIED,
  LDB,
  RDB,
  LUB,
  RUB,
  LDF,
  RDF,
  LUF,
  RUF,
}

export interface PackOptions {
  version: number;
  from: CoordinateSystem;
  sh1Bits: number;
  shRestBits: number;
  hasSafeOrbit: boolean;
  safeOrbitElevationMin: number;
  safeOrbitElevationMax: number;
  safeOrbitRadiusMin: number;
}

export interface UnpackOptions {
  to: CoordinateSystem;
}

export class GaussianCloud {
  numPoints: number;
  shDegree: number;
  antialiased: boolean;
  positions: Float32Array;
  scales: Float32Array;
  rotations: Float32Array;
  alphas: Float32Array;
  colors: Float32Array;
  sh: Float32Array;
}

// Wrappers for the module interface
export interface SpzModule {
  // Enum
  CoordinateSystem: typeof CoordinateSystem;

  // Loaders / Savers
  loadSpzFromBuffer(data: Uint8Array, options: UnpackOptions): GaussianCloud;
  saveSpzToBuffer(cloud: GaussianCloud, options: PackOptions): Uint8Array;

  LATEST_SPZ_HEADER_VERSION: number;
}

export default function createSpzModule(overrides?: any): Promise<SpzModule>;
