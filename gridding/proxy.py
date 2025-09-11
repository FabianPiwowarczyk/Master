import numpy as np

# --- helpers ---------------------------------------------------------------
def running_mean_1d(y, win=9):
    """NaN-safe centered running mean."""
    y = np.asarray(y, float)
    m = np.isfinite(y).astype(float)
    y0 = np.where(np.isfinite(y), y, 0.0)
    k = np.ones(win)
    num = np.convolve(y0, k, mode='same')
    den = np.convolve(m,  k, mode='same')
    out = num / np.where(den == 0, np.nan, den)
    return out

def vortex_weights(lat, month, inner=50, outer=60, floor=0.3):
    """
    Cosine taper to down-weight vortex latitudes in the active hemisphere.
    NH season: Dec–Apr (12,1,2,3,4);  SH season: Jun–Oct (6–10).
    """
    lat = np.asarray(lat)
    w = np.ones_like(lat, float)

    def taper(mask_sign):
        # cosine ramp between inner and outer, then clamp at 'floor'
        a = np.abs(lat)
        ramp = 0.5 * (1 + np.cos(np.pi * (a - inner) / (outer - inner)))
        wloc = np.where(a <= inner, 1.0,
               np.where(a >= outer, floor, floor + (1 - floor) * ramp))
        return np.where(mask_sign, wloc, 1.0)

    # activate NH taper in winter/spring
    if month in (12, 1, 2, 3, 4):
        w = w * taper(lat > 0)
    # activate SH taper in winter/spring
    if month in (6, 7, 8, 9, 10):
        w = w * taper(lat < 0)

    return w


def zonal_mean_safe(X):
    X = np.asarray(X, float)             # (lat, lon)
    valid = np.isfinite(X)
    cnt = valid.sum(axis=1)              # number of valid longitudes per lat
    s   = np.nansum(np.where(valid, X, 0.0), axis=1)
    zm  = np.divide(s, cnt,
                    out=np.full(cnt.shape, np.nan),
                    where=cnt > 0)       # NaN where there is no data
    return zm


def tropospheric_proxy_from_totals(X, lat, month, smooth_win=9, trop_anchor=20):
    """
    Method C: build a tropospheric proxy from a gridded total-column field X(lat,lon).
    1) remove zonal mean (zonal anomaly)
    2) down-weight vortex latitudes in season
    3) add back a *smoothed* zonal mean as a baseline
    4) anchor the tropical (|lat|<=trop_anchor) mean to the original
    """
    # 1) zonal mean and anomaly
    # assume X shape is (nlat, nlon); adjust axis if needed
    zm = zonal_mean_safe(X)                       # (nlat,)
    anom = X - zm[:, None]                             # (nlat, nlon)

    # 2) down-weight high-lat anomalies during vortex season
    wlat = vortex_weights(lat, month)                  # (nlat,)
    anom_w = anom * wlat[:, None]

    # 3) very smooth meridional baseline (keeps broad tropospheric gradient)
    zm_smooth = running_mean_1d(zm, win=smooth_win)    # (nlat,)
    proxy = anom_w + zm_smooth[:, None]

    # 4) anchor: keep the original tropical mean exactly
    trop = np.abs(lat) <= trop_anchor
    if np.any(trop):
        orig_trop_mean  = np.nanmean(X[trop, :])
        proxy_trop_mean = np.nanmean(proxy[trop, :])
        proxy += (orig_trop_mean - proxy_trop_mean)

    return proxy, zm, zm_smooth, wlat