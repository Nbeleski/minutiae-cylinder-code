#pragma once

#include "fingerprint_base.h"

#include <array>
#include <bitset>

constexpr float NS = 16.0f; // 8.0
constexpr float ND = 6.0f;
constexpr float CYLINDER_RADIUS = 70.0f;
constexpr float DELTA_S = 2.0f * CYLINDER_RADIUS / NS;
constexpr float DELTA_D = 2.0f * PI / ND;

constexpr float SIGMA_S = 28.0f / 3.0f;
constexpr float SIGMA_D = 2.0f * PI / 9.0f;
constexpr float SQRT_TWO_PI = 2.5066283f;

constexpr float OMEGA = 50.0f;

constexpr float MIN_VALID_CELLS  = 0.75f;	 // % minimo de celulas para um cilindro ser valido
constexpr int MIN_CONTRIBUTING_MINUTIAE = 2; // numero minimo de minucias contribuintes para um cilindro valido
constexpr float MAX_ANG_DIFF = PI / 2.0f;	 // angulo maximo na orientação da minucia/cilindro

constexpr float NP_MIN = 4.0f;
constexpr float NP_MAX = 12.0f;

constexpr size_t CYLINDER_SIZE = (NS * NS * ND);

// bitwise params
constexpr float MICRO_PSI = 0.01f;
typedef std::bitset<CYLINDER_SIZE> LinearizedBitCylinder;

