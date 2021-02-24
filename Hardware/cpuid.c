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

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/sysinfo.h>

#include "cpuid.h"
#include "libmsr_error.h"

void cpuid(uint64_t leaf, uint64_t *rax, uint64_t *rbx, uint64_t *rcx, uint64_t *rdx)
{
    asm volatile (
        "\txchg %%rbx, %%rdi\n"
        "\tcpuid\n"
        "\txchg %%rbx, %%rdi"
        : "=a" (*rax), "=D" (*rbx), "=c" (*rcx), "=d" (*rdx)
        : "a" (leaf)
    );
}

void cpuid_detect_core_conf(uint64_t *coresPerSocket, uint64_t *hyperThreads, uint64_t *sockets, int *HTenabled)
{
    unsigned first = 0;
    unsigned second = 0;
    uint64_t rax = 0xb;
    uint64_t rbx = 0;
    uint64_t rcx = 0x0;
    uint64_t rdx = 0;
    static int init = 0;
    int allcores = 0;

    // Use rcx = 0 to see if hyperthreading is supported. If > 1, then there is
    // HT.
    // Use rcx = 1 to see how many cores are available per socket (including
    // HT, if supported).
    if (!init)
    {
        FILE *thread;
        thread = fopen("/sys/devices/system/cpu/cpu0/topology/thread_siblings_list", "r");
        if (thread == NULL)
        {
            libmsr_error_handler("cpuid_detect_core_conf(): Unable to open /sys/devices/system/cpu/cpu0/topology/thread_siblings_list", LIBMSR_ERROR_PLATFORM_ENV, getenv("HOSTNAME"), __FILE__, __LINE__);
        }
        int ret = fscanf(thread, "%u,%u", &first, &second);
        if (ret < 2 || ret == EOF)
        {
            /* Hyperthreading is disabled. */
            *HTenabled = 0;
        }
        else
        {
            *HTenabled = 1;
        }
        if (thread != NULL)
        {
            fclose(thread);
        }
        init = 1;
    }

    asm volatile(
        "cpuid"
        : "=a" (rax), "=b" (rbx), "=c" (rcx), "=d" (rdx)
        : "0" (rax), "2"(rcx)
    );
    *hyperThreads = ((rbx) & 0xFFFF);
    rax = 0xb;
    rbx = 0;
    rcx = 0x1;
    rdx = 0;

    asm volatile(
        "cpuid"
        : "=a" (rax), "=b" (rbx), "=c" (rcx), "=d" (rdx)
        : "0" (rax), "2"(rcx)
    );
    *coresPerSocket = ((rbx) & 0xFFFF) / *hyperThreads;
    // get_nprocs_conf() returns max number of logical processors (including
    // hyperthreading)
    // get_nprocs() returns num logical processors depending on whether
    // hyperthreading is enabled or not
    allcores = get_nprocs_conf();
    if (allcores == *coresPerSocket * (*HTenabled + 1))
    {
        *sockets = 1;
    }
    else
    {
        *sockets = 2;
    }
#ifdef CPUID_DEBUG
    fprintf(stderr, "%s::%d DEBUG: allcores is %d, and register has 0x%lx\n", __FILE__, __LINE__, allcores, rbx);
    fprintf(stderr, "%s::%d DEBUG: hyper threads is %ld, cores per socket is %ld, sockets is %ld, HT is %d\n", __FILE__, __LINE__, *hyperThreads, *coresPerSocket, *sockets, *HTenabled);
#endif
}

void cpuid_get_model(uint64_t *model)
{
    /* Set rax to 1 which indicates we want processor info and feature bits. */
    uint64_t rax = 1;
    uint64_t rbx = 0;
    uint64_t rcx = 0;
    uint64_t rdx = 0;

    /* This is how the linux kernel does it. */
    asm volatile (
        "cpuid"
        : "=a" (rax), "=b" (rbx), "=c" (rcx), "=d" (rdx)
        : "0" (rax), "2"(rcx)
    );

#ifdef CPUID_DEBUG
    fprintf(stderr, "%s::%d DEBUG: rax is %lx\n", __FILE__, __LINE__, rax);
    fprintf(stderr, "%s::%d DEBUG: model is %lx, family is %lx, extd model is %lx, extd family ix %lx\n", __FILE__, __LINE__, (rax >> 4) & 0xF, (rax >> 8) & 0xF, (rax >> 16) & 0xF, (rax >> 20) & 0xFF);
#endif
    *model = ((rax >> 4) & 0xF) | ((rax >> 12) & 0xF0);
}

void cpuidInput_rax_rcx(uint64_t leafa, uint64_t leafc, uint64_t *rax, uint64_t *rbx, uint64_t *rcx, uint64_t *rdx)
{
    asm volatile(
        "xchg %%rbx, %%rdi\n"
        "\tcpuid\n"
        "\txchg %%rbx, %%rdi"
        : "=a" (*rax), "=D" (*rbx), "=c" (*rcx), "=d" (*rdx)
        : "a" (leafa), "c" (leafc)
    );
}

bool cpuid_mperf_and_aperf_avail(void)
{
    /* See Manual Vol 3B, Section 14.2 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rcx, 0, 0) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_timestamp_counter_avail(void)
{
    /* See Manual Vol 3B, Section 17.15 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 1;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rdx, 4, 4) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

// NOTE: DS_AREA: does not need to check cpuid specifically (Manual Vol 3C) may
// be in another source
// NOTE: PEBS_ENABLE: does not need to check cpuid specifically (Manual Vol 3C)
// may be in another source

int cpuid_num_pmc(void)
{
    /* See Manual Vol 3B, Section 18.2.5 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 10; // 0A

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    // If CPUID.0AH:rax[15:8] > 3, then up to PMC3 is usable
    // If CPUID.0AH:rax[15:8] > 2, then up to PMC2 is usable
    // If CPUID.0AH:rax[15:8] > 1, then up to PMC1 is usable
    // If CPUID.0AH:rax[15:8] > 0, then only PMC0 is usable
    // If CPUID.0AH:rax[15:8] == 0, then none are usable
    return MASK_VAL(rax, 15, 8);
}

int cpuid_num_perfevtsel(void)
{
    /* See Manual Vol 3B, Section 18.2.1.1 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 10; // 0A

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    // If CPUID.0AH:rax[15:8] > 3, then up to PerfEvtSel3 is usable
    // If CPUID.0AH:rax[15:8] > 2, then up to PerfEvtSel2 is usable
    // If CPUID.0AH:rax[15:8] > 1, then up to PerfEvtSel1 is usable
    // If CPUID.0AH:rax[15:8] > 0, then only PerfEvtSel0 is usable
    // If CPUID.0AH:rax[15:8] == 0, then none are usable
    return MASK_VAL(rax, 15, 8);
}

bool cpuid_perf_global_ctrl_EN_PMC(void)
{
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 10; // 0A

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 7, 0) > 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_perf_global_ctrl_EN_FIXED_CTRnum(void)
{
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 10; // 0A

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 7, 0) > 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_misc_enable_TurboModeDisable(void)
{
    /* See Manual Vol 3C, Table 35-12 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 1, 1) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_misc_enable_xTPRMessageDisable(void)
{
    /* See Manual Vol 3C, Table 35-2 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 1;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rcx, 14, 14) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_misc_enable_XDBitDisable(void)
{
    /* See Manual Vol 3C, Table 35-2 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 2147483649; // 80000001H

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rdx, 20, 20) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_enable_ExtendedClockMod(void)
{
    /* See Manual Vol 3B, Section 14.7.3.1 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 5, 5) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_therm_status_enable_ThermalThresholds(void)
{
    /* See Manual Vol 3B, Section 14.7.2.2 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 1;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rcx, 8, 8) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_therm_stat_enable_PowerLimitNotify(void)
{
    /* See Manual Vol 3B, Section 14.7.6 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 4, 4) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_therm_status_enable_DigitalReadout(void)
{
    /* See Manual Vol 3B, Section 14.7.5.1 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 0, 0) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_therm_interrupt_enable_PowerLimitNotify(void)
{
    /* See Manual Vol 3B, Section 14.7.6 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 4, 4) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool cpuid_pkg_therm_enable_status_and_interrupt(void)
{
    /* See Manual Vol 3B, Section 14.8 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 6;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    if (MASK_VAL(rax, 6, 6) == 1)
    {
        return true;
    }
    else
    {
        return false;
    }
}

uint64_t cpuid_MaxLeaf(void)
{
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 0;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    return rax;
}

void cpuid_printVendorID(void)
{
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 0;
    int i = 0;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    for (i = 0; i < 32; i += 8)
    {
        printf("%c", (int)MASK_VAL(rbx, 7+i, 0+i));
    }
    for (i = 0; i < 32; i += 8)
    {
        printf("%lc", (int)MASK_VAL(rdx, 7+i, 0+i));
    }
    for (i = 0; i < 32; i += 8)
    {
        printf("%lc", (int)MASK_VAL(rcx, 7+i, 0+i));
    }
    printf("\n");
}

int cpuid_pkg_MaxPhysicalProcessorCores(void)
{
    /* See Manual Vol 3A, Section 8.6 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leafa = 4;
    int leafc = 0;

    cpuidInput_rax_rcx(leafa, leafc, &rax, &rbx, &rcx, &rdx);
    return MASK_VAL(rax, 31, 26) + 1;
}

int cpuid_pkg_MaxLogicalProcessors(void)
{
    /* See Manual Vol 3A, Section 8.6 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 1;

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    return MASK_VAL(rbx, 23, 16);
}

int cpuid_num_fixed_counters(void)
{
    /* See Manual Vol 3B, Section 18.2.2.1 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 10; //0A

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    return MASK_VAL(rdx, 4, 0);
}

int cpuid_width_fixed_counters(void)
{
    /* See Manual Vol 3B, Section 18.2.2.1 for details. */
    uint64_t rax, rbx, rcx, rdx;
    int leaf = 10; //0A

    cpuid(leaf, &rax, &rbx, &rcx, &rdx);
    return MASK_VAL(rdx, 12, 5);
}
