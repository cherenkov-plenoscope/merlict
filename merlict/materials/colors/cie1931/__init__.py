from . import rgb_matching_functions
from . import standard_observer

import numpy as np


def visible_wavelength_range():
    return [400e-9, 700e-9]


def visible_wavelengths(num_steps_in_visible_range=31):
    w = [200e-9]
    vstart, vstop = visible_wavelength_range()
    visible = np.linspace(vstart, vstop, num_steps_in_visible_range)
    w += visible.tolist()
    w += [1200e-9]
    return w


class Cie1931:
    """
    The CIE 1931 colorspace.

    The geometry in merlict desribes its material properties using spectral
    functions vs. the wavelength of light.

    This tool helps you to create spectra to match a specific color perception.
    But watch out, the spectra created here are only of visual value and might
    not represent any existing materrial.
    """

    def __init__(self, wavelengths=None):
        if wavelengths is None:
            self.wavelengths = visible_wavelengths()
        else:
            self.wavelengths = wavelengths

        o_x = standard_observer.x()
        o_y = standard_observer.y()
        o_z = standard_observer.z()
        self.x = np.interp(x=self.wavelengths, xp=o_x[:, 0], fp=o_x[:, 1])
        self.y = np.interp(x=self.wavelengths, xp=o_y[:, 0], fp=o_y[:, 1])
        self.z = np.interp(x=self.wavelengths, xp=o_z[:, 0], fp=o_z[:, 1])
        self.x_sum = np.sum(self.x)
        self.y_sum = np.sum(self.y)
        self.z_sum = np.sum(self.z)

        m_r = rgb_matching_functions.r()
        m_g = rgb_matching_functions.g()
        m_b = rgb_matching_functions.b()
        self.r = np.interp(x=self.wavelengths, xp=m_r[:, 0], fp=m_r[:, 1])
        self.g = np.interp(x=self.wavelengths, xp=m_g[:, 0], fp=m_g[:, 1])
        self.b = np.interp(x=self.wavelengths, xp=m_b[:, 0], fp=m_b[:, 1])
        self.r_sum = np.sum(self.r)
        self.g_sum = np.sum(self.g)
        self.b_sum = np.sum(self.b)

        wvl_start, wvl_stop = visible_wavelength_range()
        self.visible_wavelength_start = wvl_start
        self.visible_wavelength_stop = wvl_stop

    def wavelength_is_visible(self, wavelength):
        return (
            wavelength > self.visible_wavelength_start
            and wavelength < self.visible_wavelength_stop
        )

    def _white_spectrum_visible_ones_invisible_zeros(self):
        wvl_start, wvl_stop = visible_wavelength_range()
        spectrum = np.ones(len(self.wavelengths))
        for w in range(len(self.wavelengths)):
            wvl = self.wavelengths[w]
            if not self.wavelength_is_visible(wvl):
                spectrum[w] = 0.0
        return spectrum

    def approximate_spectrum_xyz(self, xyz, max_num_iterations=100, mag=0.05):
        """
        This tries to make the spectrum as borad as possible to avoid
        monochromatic lines.

        Monochromatic lines are the easy way to get the spectrum we want but it
        is not ideal if the light source has not a broad spectrum itself.

        Parameters
        ----------
        xyz : arraylike floats
            The xyz response of the 'standard observer'.

        Returns
        -------
        spectrum : arraylike floats shape (N,2)
            wavelength/m in [:, 0], amplitudes in [:, 1].
            N pairs of wavelength and amplitude.
        """
        _spectrum = self._approximate_spectrum(
            val=xyz, mode="xyz", max_num_iterations=max_num_iterations, mag=mag
        )
        return np.asarray([self.wavelengths, _spectrum]).T

    def approximate_spectrum_rgb(self, rgb, max_num_iterations=100, mag=0.05):
        _spectrum = self._approximate_spectrum(
            val=rgb, mode="rgb", max_num_iterations=max_num_iterations, mag=mag
        )
        return np.asarray([self.wavelengths, _spectrum]).T

    def _approximate_spectrum(
        self, val, mode="xyz", max_num_iterations=100, mag=0.05
    ):
        """
        A crude iterative search of a spectrum that best yields the requested
        response 'xyz' of the standard observer.
        """
        if mode == "xyz":
            observe_spectrum = self._observe_spectrum_xyz
        elif mode == "rgb":
            observe_spectrum = self._observe_spectrum_rgb
        else:
            assert False, "Expected mode 'xyz' or 'rgb'."

        val /= np.mean(val)

        spectrum = self._white_spectrum_visible_ones_invisible_zeros()

        wvl_start, wvl_stop = visible_wavelength_range()

        last_delta = 2 * np.sum(val)

        for i in range(max_num_iterations):
            _val = observe_spectrum(spectrum)
            delta = np.linalg.norm(val - _val)
            if last_delta < delta * (1 + 1e-3):
                break
            last_delta = delta

            for w in range(len(self.wavelengths)):
                wvl = self.wavelengths[w]
                if not self.wavelength_is_visible(wvl):
                    continue

                _spectrum = spectrum.copy()

                _deltas = [delta]
                _signs = [0]

                _spectrum[w] = (1 - mag) * spectrum[w]
                _val = observe_spectrum(_spectrum)
                _delta = np.linalg.norm(val - _val)
                _deltas.append(_delta)
                _signs.append(-1)

                _spectrum[w] = (1 + mag) * spectrum[w]
                _val = observe_spectrum(_spectrum)
                _delta = np.linalg.norm(val - _val)
                _deltas.append(_delta)
                _signs.append(+1)

                if np.std(_deltas) < 1e-6:
                    best_sign = 0
                else:
                    best_case = np.argmin(_deltas)
                    best_sign = _signs[best_case]

                if best_sign == 0:
                    pass
                elif best_sign == -1:
                    spectrum[w] = spectrum[w] * (1 - mag)
                elif best_sign == +1:
                    spectrum[w] = spectrum[w] * (1 + mag)
        return spectrum / spectrum.max()

    def _observe_spectrum_xyz(self, spectrum):
        xi = np.sum(spectrum * self.x) / self.x_sum
        yi = np.sum(spectrum * self.y) / self.y_sum
        zi = np.sum(spectrum * self.z) / self.z_sum
        return np.asarray([xi, yi, zi])

    def observe_spectrum_xyz(self, spectrum):
        """
        Parameters
        ----------
        spectrum : arraylike floats shape (N,2)
            wavelength/m in [:, 0], amplitudes in [:, 1].
            N pairs of wavelength and amplitude.

        Returns
        -------
        xyz : arraylike floats
            The xyz response of the 'standard observer'.
        """
        _spectrum = np.interp(
            x=self.wavelengths, xp=spectrum[:, 0], fp=spectrum[:, 1]
        )
        return self._observe_spectrum_xyz(_spectrum)

    def _observe_spectrum_rgb(self, spectrum):
        r = np.sum(spectrum * self.r) / self.r_sum
        g = np.sum(spectrum * self.g) / self.g_sum
        b = np.sum(spectrum * self.b) / self.b_sum
        return np.asarray([r, g, b])

    def observe_spectrum_rgb(self, spectrum):
        _spectrum = np.interp(
            x=self.wavelengths, xp=spectrum[:, 0], fp=spectrum[:, 1]
        )
        return self._observe_spectrum_rgb(_spectrum)

    def _print_normalized(self, channel):
        ox = self.x.copy()
        oy = self.y.copy()
        oz = self.z.copy()
        m = np.max([np.max(ox), np.max(oy), np.max(oz)])
        ox /= m
        oy /= m
        oz /= m

        s = ""
        for cc in [("x", ox), ("y", oy), ("z", oz)]:
            s += f"{cc[0]:s} = {{\n"
            for w in range(len(self.wavelengths)):
                wvl = self.wavelengths[w]
                val = cc[1][w]
                s += f"    {{{wvl:.3e}, {val:.3e}}},\n"
            s += "}\n"
        return s
