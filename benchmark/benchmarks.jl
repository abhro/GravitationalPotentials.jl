# SPDX-License-Identifier: MIT

using GravitationalPotentials
using BenchmarkTools

SUITE = BenchmarkGroup()
SUITE["rand"] = @benchmarkable rand(10)

# Write your benchmarks here.
