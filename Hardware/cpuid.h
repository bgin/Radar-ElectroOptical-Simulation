/*
 * Copyright (c) 2013-2017, Lawrence Livermore National Security, LLC.
 *
 * Produced at the Lawrence Livermore National Laboratory. Written by:
 *     Barry Rountree <rountree@llnl.gov>,
 *     Scott Walker <walker91@llnl.gov>, and
 *     Kathleen Shoga <shoga1@llnl.gov>.
 *
 * LLNL-CODE-645430
 *
 * All rights reserved.
 *
 * This file is part of libmsr. For details, see https://github.com/LLNL/libmsr.git.
 *
 * Please also read libmsr/LICENSE for our notice and the LGPL.
 *
 * libmsr is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License (as published by the Free
 * Software Foundation) version 2.1 dated February 1999.
 *
 * libmsr is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU Lesser
 * General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with libmsr; if not, write to the Free Software Foundation, Inc., 59
 * Temple Place, Suite 330, Boston, MA 02111-1307 USA.
 */

#ifndef CPUID_H_INCLUDE
#define CPUID_H_INCLUDE

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

/// @brief Create a mask from bit m to n (63 >= m >= n >= 0).
///
/// Example: MASK_RANGE(4,2) --> (((1<<((4)-(2)+1))-1)<<(2))
///                              (((1<<          3)-1)<<(2))
///                              ((               4-1)<<(2))
///                              (                  3)<<(2))
///                              (                       24) = b11000
#define MASK_RANGE(m,n) ((((uint64_t)1<<((m)-(n)+1))-1)<<(n))

/// @brief Return the value of x after applying bitmask (m,n) (63 >= m >= n >= 0).
///
/// Example: MASK_RANGE(17,4,2) --> 17&24 = b10001 & b11000 = b10000
#define MASK_VAL(x,m,n) (((uint64_t)(x)&MASK_RANGE((m),(n)))>>(n))

/// @brief Get processor identification and feature information to the 32-bit
/// EAX, EBX, ECX, and EDX registers.
///
/// @param [in] leaf Level of basic processor information to be returned.
///
/// @param [out] rax EAX register contents.
///
/// @param [out] rbx EBX register contents.
///
/// @param [out] rcx ECX register contents.
///
/// @param [out] rdx EDX register contents.
void cpuid(uint64_t leaf,
           uint64_t *rax,
           uint64_t *rbx,
           uint64_t *rcx,
           uint64_t *rdx);

/// @brief Determine platform configuration (i.e. number of sockets, number of
/// cores/socket, number of logical processors, etc.).
///
/// @param [out] coresPerSocket Number of cores per socket (including
/// hyperthreads, if enabled).
///
/// @param [out] hyperThreads Number of threads including hyperthreads.
///
/// @param [out] sockets Number of sockets.
///
/// @param [out] HTenabled Hyperthreading enabled/disabled.
void cpuid_detect_core_conf(uint64_t *coresPerSocket,
                            uint64_t *hyperThreads,
                            uint64_t *sockets,
                            int *HTenabled);

/// @brief Retrieve the CPU model number.
///
/// @param [out] model CPU model number.
void cpuid_get_model(uint64_t *model);

/// @brief Get max number of physical processor cores on the platform.
///
/// @param [in] leafa Initial value for EAX.
///
/// @param [in] leafc Initial value for ECX.
///
/// @param [out] rax EAX register contents.
///
/// @param [out] rbx EBX register contents.
///
/// @param [out] rcx ECX register contents.
///
/// @param [out] rdx EDX register contents.
void cpuidInput_rax_rcx(uint64_t leafa,
                        uint64_t leafc,
                        uint64_t *rax,
                        uint64_t *rbx,
                        uint64_t *rcx,
                        uint64_t *rdx);

/// @brief Check availability of p-state hardware coordination feedback,
/// indicating presence of IA32_APERF and IA32_MPERF.
///
/// @return True if IA32_APERF and IA32_MPERF are available, else false.
bool cpuid_mperf_and_aperf_avail(void);

/// @brief Check availability of the time-stamp counter.
///
/// @return True if the time-stamp counter is available, else false.
bool cpuid_timestamp_counter_avail(void);

/*****************************************/
/* Performance Monitoring Counters (PMC) */
/*****************************************/

/// @brief Determine which performance monitoring counters (PMCs) are available.
///
/// If greater than 3, then up to PMC3 is usable. If greater than 2, then up to
/// PMC2 is usable. If greater than 1, then up to PMC1 is usable. If greater than
/// 0, then only PMC0 is usable. If equal to 0, then no PMCs are usable.
///
/// @return Number of PMCs are available.
int cpuid_num_pmc(void);

/*****************************************/
/* Performance Event Select (PerfEvtSel) */
/* (0x186, 0x187, 0x188, 0x189)          */
/*****************************************/

/// @brief Determine which performance event select MSRs (PerfEvtSelX) are
/// available.
///
/// If greater than 3, then up to PerfEvtSel3 is usable. If greater than 2,
/// then up to PerfEvtSel2 is usable. If greater than 1, then up to PerfEvtSel1
/// is usable. If greater than 0, then only PerfEvtSel0 is usable. If equal to
/// 0, then no PerfEvtSel MSRs are usable.
///
/// @return Number of PerfEvtSel MSRs are available.
int cpuid_num_perfevtsel(void);

/***************************/
/* PERF_GLOBAL_CTL (0x38f) */
/***************************/

/// @brief Determine if performance monitoring capabilities are available.
///
/// @return True if architectural performance monitoring capabilities are
/// supported, else false.
bool cpuid_perf_global_ctrl_EN_PMC(void);

/// @brief Check availability of the three fixed-function performance counters:
/// Instr_Retired.Any, CPU_CLK_Unhalted.Core, and CPU_CLK_Unhalted.Ref.
///
/// @return True if IA32_FIXED_CTR0, IA32_FIXED_CTR1, and IA32_FIXED_CTR2
/// exist, else false.
bool cpuid_perf_global_ctrl_EN_FIXED_CTRnum(void);

/***************/
/* MISC_ENABLE */
/***************/

/// @brief Check if Turbo mode is disabled
///
/// @return True if Turbo mode is disabled, else false.
bool cpuid_misc_enable_TurboModeDisable(void);

/// @brief Check if xTPR messages are disabled.
///
/// @return True if xTPR messages are disabled, else false.
bool cpuid_misc_enable_xTPRMessageDisable(void);

/// @brief Check if XD bit is disabled.
///
/// @return True if XD bit is disabled, else false.
bool cpuid_misc_enable_XDbitDisable(void);

/*************/
/* CLOCK_MOD */
/*************/

/// @brief Check if extended on-demand clock modulation is enabled.
///
/// @return True if extended on-demand clock modulation is enabled, else false.
bool cpuid_enable_ExtendedClockMod(void);

/*********************/
/* THERMAL Functions */
/*********************/

/// @brief Check if Thermal Threshold #1 status and log and Thermal Threshold
/// #2 status and log are enabled in IA32_THERM_STATUS.
///
/// @return True if TT#1 and TT#2 status and log are enabled, else false.
bool cpuid_therm_status_enable_ThermalThresholds(void);

/// @brief Check if power limitation status and log are enabled in
/// IA32_THERM_STATUS.
///
/// @return True if power limitation status and log are enabled, else false.
bool cpuid_therm_status_enable_PowerLimitNotify(void);

/// @brief Check if digital readout is enabled in IA32_THERM_STATUS.
///
/// @return True if digital readout is enabled, else false.
bool cpuid_therm_status_enable_DigitalReadout(void);

/// @brief Check if power limit notification enable is enabled in
/// IA32_THERM_INTERRUPT.
///
/// @return True if power limitation status and log are enabled, else false.
bool cpuid_therm_interrupt_enable_PowerLimitNotify(void);

/// @brief Check if IA32_PACKAGE_THERM_STATUS and IA32_PACKAGE_THERM_INTERRUPT
/// are enabled on the platform.
///
/// @return True if package-level thermal status and interrupt register are
/// enabled, else false.
bool cpuid_pkg_therm_enable_status_and_interrupt(void);

/************************/
/* General Machine Info */
/************************/

/// @brief Get highest EAX value (i.e., leaf) supported by the platform.
///
/// @return Highest EAX leaf supported by the platform.
uint64_t cpuid_MaxLeaf(void);

/// @brief Print platform vendor ID.
void cpuid_printVendorID(void);

/// @brief Get max number of addressable IDs attributable to processor cores
/// in the physical package.
///
/// @return Max number of processor cores in the package.
int cpuid_pkg_MaxPhysicalProcessorCores(void);

/// @brief Get max number of addressable IDs for logical processors in a
/// physical package.
///
/// @return Max number of logical processors in the package.
int cpuid_pkg_MaxLogicalProcessors(void);

/// @brief Get number of fixed-function performance counters on the platform.
///
/// @return Number of available fixed-function performance counters.
int cpuid_num_fixed_counters(void);

/// @brief Get bit width of fixed-function performance counters on the platform.
///
/// @return Bit width of fixed-function performance counters.
int cpuid_width_fixed_counters(void);

#ifdef __cplusplus
}
#endif
#endif
