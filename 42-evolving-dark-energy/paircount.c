/*
 * Weighted pair counts in bins of separation s and |mu| (the cosine of the
 * angle between the pair and the line of sight), for the two-point
 * correlation function. Written for this project; compiled on first use by
 * bao_measurement.py (gcc -O3 -fopenmp -shared).
 *
 * Points are sorted into cubic cells of side >= s_max / 2, so every pair
 * closer than s_max lies within two cells of each other in each direction.
 * Each pair is counted once (auto) or once per (i in A, j in B) (cross).
 *
 * Besides the total, the counts of pairs touching each jackknife region are
 * kept: leaving region k out removes exactly those pairs, so the delete-one
 * jackknife counts are total - touch[k], with no recounting.
 *
 * Also accumulated: the weighted mean redshift of the pairs in each s bin
 * (sum of w_i w_j (v_i + v_j) / 2), from which z_eff is computed.
 */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>

typedef struct {
    int nx, ny, nz;
    double x0, y0, z0, size;
    long *start;   /* first point of each cell in the sorted order (ncell + 1) */
    long *order;   /* point indices sorted by cell */
} Grid;

static long cell_of(const Grid *g, double x, double y, double z) {
    int i = (int)((x - g->x0) / g->size), j = (int)((y - g->y0) / g->size),
        k = (int)((z - g->z0) / g->size);
    if (i >= g->nx) i = g->nx - 1;
    if (j >= g->ny) j = g->ny - 1;
    if (k >= g->nz) k = g->nz - 1;
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (k < 0) k = 0;
    return ((long)i * g->ny + j) * g->nz + k;
}

static void grid_build(Grid *g, const double *pos, long n) {
    long ncell = (long)g->nx * g->ny * g->nz;
    long *count = calloc(ncell + 1, sizeof(long));
    long *cell = malloc(n * sizeof(long));
    for (long p = 0; p < n; p++) {
        cell[p] = cell_of(g, pos[3 * p], pos[3 * p + 1], pos[3 * p + 2]);
        count[cell[p] + 1]++;
    }
    g->start = malloc((ncell + 1) * sizeof(long));
    g->start[0] = 0;
    for (long c = 0; c < ncell; c++) g->start[c + 1] = g->start[c] + count[c + 1];
    long *fill = calloc(ncell, sizeof(long));
    g->order = malloc(n * sizeof(long));
    for (long p = 0; p < n; p++) g->order[g->start[cell[p]] + fill[cell[p]]++] = p;
    free(count);
    free(cell);
    free(fill);
}

/*
 * pos1 (n1 x 3), w1, reg1, v1: first catalogue; pos2 ... : second (ignored if auto).
 * Output arrays (zeroed by the caller):
 *   total[ns * nmu], touch[nreg * ns * nmu], vsum[ns], wsum[ns]
 */
int count_pairs(const double *pos1, const double *w1, const int *reg1, const double *v1, long n1,
                const double *pos2, const double *w2, const int *reg2, const double *v2, long n2,
                int autocorr, double s_max, int ns, int nmu, int nreg,
                double *total, double *touch, double *vsum, double *wsum, int nthreads) {
    if (autocorr) { pos2 = pos1; w2 = w1; reg2 = reg1; v2 = v1; n2 = n1; }
    /* bounding box of both catalogues */
    double lo[3] = {1e300, 1e300, 1e300}, hi[3] = {-1e300, -1e300, -1e300};
    for (long p = 0; p < n1; p++)
        for (int a = 0; a < 3; a++) {
            if (pos1[3 * p + a] < lo[a]) lo[a] = pos1[3 * p + a];
            if (pos1[3 * p + a] > hi[a]) hi[a] = pos1[3 * p + a];
        }
    for (long p = 0; p < n2; p++)
        for (int a = 0; a < 3; a++) {
            if (pos2[3 * p + a] < lo[a]) lo[a] = pos2[3 * p + a];
            if (pos2[3 * p + a] > hi[a]) hi[a] = pos2[3 * p + a];
        }
    Grid g;
    g.size = s_max / 2.0;
    g.x0 = lo[0]; g.y0 = lo[1]; g.z0 = lo[2];
    g.nx = (int)((hi[0] - lo[0]) / g.size) + 1;
    g.ny = (int)((hi[1] - lo[1]) / g.size) + 1;
    g.nz = (int)((hi[2] - lo[2]) / g.size) + 1;
    Grid g2 = g;
    grid_build(&g, pos1, n1);
    if (autocorr) g2 = g; else grid_build(&g2, pos2, n2);

    const double s2max = s_max * s_max, ds = s_max / ns;
    const long nbin = (long)ns * nmu;
    const long ncell = (long)g.nx * g.ny * g.nz;
    if (nthreads > 0) omp_set_num_threads(nthreads);

#pragma omp parallel
    {
        double *t_total = calloc(nbin, sizeof(double));
        double *t_touch = calloc((long)nreg * nbin, sizeof(double));
        double *t_vsum = calloc(ns, sizeof(double));
        double *t_wsum = calloc(ns, sizeof(double));
#pragma omp for schedule(dynamic, 16)
        for (long c1 = 0; c1 < ncell; c1++) {
            long a0 = g.start[c1], a1 = g.start[c1 + 1];
            if (a0 == a1) continue;
            int i1 = (int)(c1 / ((long)g.ny * g.nz)), j1 = (int)((c1 / g.nz) % g.ny),
                k1 = (int)(c1 % g.nz);
            for (int di = -2; di <= 2; di++) {
                int i2 = i1 + di;
                if (i2 < 0 || i2 >= g.nx) continue;
                for (int dj = -2; dj <= 2; dj++) {
                    int j2 = j1 + dj;
                    if (j2 < 0 || j2 >= g.ny) continue;
                    for (int dk = -2; dk <= 2; dk++) {
                        int k2 = k1 + dk;
                        if (k2 < 0 || k2 >= g.nz) continue;
                        /* cells whose closest points are farther than s_max: skip */
                        double gx = abs(di) > 1 ? (abs(di) - 1) * g.size : 0,
                               gy = abs(dj) > 1 ? (abs(dj) - 1) * g.size : 0,
                               gz = abs(dk) > 1 ? (abs(dk) - 1) * g.size : 0;
                        if (gx * gx + gy * gy + gz * gz > s2max) continue;
                        long c2 = ((long)i2 * g.ny + j2) * g.nz + k2;
                        if (autocorr && c2 < c1) continue;
                        long b0 = g2.start[c2], b1 = g2.start[c2 + 1];
                        for (long a = a0; a < a1; a++) {
                            long p = g.order[a];
                            double x1 = pos1[3 * p], y1 = pos1[3 * p + 1], z1 = pos1[3 * p + 2];
                            double wp = w1[p], vp = v1[p];
                            int rp = reg1[p];
                            long bstart = (autocorr && c2 == c1) ? a + 1 : b0;
                            for (long b = bstart; b < b1; b++) {
                                long q = g2.order[b];
                                double dx = pos2[3 * q] - x1, dy = pos2[3 * q + 1] - y1,
                                       dz = pos2[3 * q + 2] - z1;
                                double s2 = dx * dx + dy * dy + dz * dz;
                                if (s2 >= s2max || s2 <= 0.0) continue;
                                double s = sqrt(s2);
                                double lx = x1 + pos2[3 * q], ly = y1 + pos2[3 * q + 1],
                                       lz = z1 + pos2[3 * q + 2];
                                double mu = fabs(dx * lx + dy * ly + dz * lz) /
                                            (s * sqrt(lx * lx + ly * ly + lz * lz));
                                int is = (int)(s / ds), im = (int)(mu * nmu);
                                if (is >= ns) is = ns - 1;
                                if (im >= nmu) im = nmu - 1;
                                double ww = wp * w2[q];
                                long bin = (long)is * nmu + im;
                                t_total[bin] += ww;
                                t_touch[(long)rp * nbin + bin] += ww;
                                if (reg2[q] != rp) t_touch[(long)reg2[q] * nbin + bin] += ww;
                                t_vsum[is] += ww * 0.5 * (vp + v2[q]);
                                t_wsum[is] += ww;
                            }
                        }
                    }
                }
            }
        }
#pragma omp critical
        {
            for (long b = 0; b < nbin; b++) total[b] += t_total[b];
            for (long b = 0; b < (long)nreg * nbin; b++) touch[b] += t_touch[b];
            for (int b = 0; b < ns; b++) { vsum[b] += t_vsum[b]; wsum[b] += t_wsum[b]; }
        }
        free(t_total); free(t_touch); free(t_vsum); free(t_wsum);
    }
    free(g.start); free(g.order);
    if (!autocorr) { free(g2.start); free(g2.order); }
    return 0;
}
