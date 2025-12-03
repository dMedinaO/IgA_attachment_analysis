from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, brentq
from scipy.stats import norm
from typing import Tuple

class AuxiliarFunctions:

    @classmethod
    def four_param_logistic(
        cls, x: np.ndarray | float, a: float, b: float, c: float, d: float
    ) -> np.ndarray:
        """
        Four-parameter logistic (4PL) model for ELISA curves.

        The function is defined as:

            y = d + (a - d) / (1 + (x / c)**b)

        where
        -----
        a : float
            Upper asymptote (maximum OD).
        b : float
            Slope (steepness) parameter.
        c : float
            Inflection point (concentration at half-maximal response).
        d : float
            Lower asymptote (minimum OD).

        Parameters
        ----------
        x : array-like or float
            Concentration values.
        a, b, c, d : float
            Model parameters.

        Returns
        -------
        numpy.ndarray
            Predicted OD values for the given concentrations.
        """
        x = np.asarray(x, dtype=float)
        return d + (a - d) / (1.0 + (x / c) ** b)
    
    @classmethod
    def prediction_band(
        cls,
        x: np.ndarray,
        params: np.ndarray,
        cov: np.ndarray,
        alpha: float = 0.05,
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Compute the fitted 4PL curve and its prediction band.

        This uses a first-order error propagation based on the Jacobian of the
        4PL model with respect to its parameters and the covariance matrix
        returned by `scipy.optimize.curve_fit`.

        Parameters
        ----------
        x : numpy.ndarray
            Concentration values at which to evaluate the curve.
        params : numpy.ndarray
            Fitted parameters [a, b, c, d].
        cov : numpy.ndarray
            Parameter covariance matrix from `curve_fit`.
        alpha : float, optional
            Significance level (default 0.05 for 95% band).

        Returns
        -------
        y : numpy.ndarray
            Predicted OD values at `x`.
        lower : numpy.ndarray
            Lower bound of the prediction band.
        upper : numpy.ndarray
            Upper bound of the prediction band.
        """
        x = np.asarray(x, dtype=float)
        y = cls.four_param_logistic(x, *params)

        # Numerical Jacobian: dy / d(theta_i)
        eps = np.sqrt(np.finfo(float).eps)
        n_params = len(params)
        J = np.empty((x.size, n_params))

        for i in range(n_params):
            dp = np.zeros_like(params)
            h = eps * max(1.0, abs(params[i]))
            dp[i] = h
            y_plus = cls.four_param_logistic(x, *(params + dp))
            y_minus = cls.four_param_logistic(x, *(params - dp))
            J[:, i] = (y_plus - y_minus) / (2.0 * h)

        # Variance of prediction: J * cov * J^T
        var_y = np.einsum("ij,jk,ik->i", J, cov, J)
        sigma_y = np.sqrt(var_y)

        z = norm.ppf(1.0 - alpha / 2.0)
        lower = y - z * sigma_y
        upper = y + z * sigma_y

        return y, lower, upper
    
    @classmethod
    def invert_4pl(
        cls,
        y_obs: float,
        params: np.ndarray,
        x_min: float = 0.0,
        x_max: float = 600.0,
    ) -> float:
        """
        Estimate the concentration for a given OD using the 4PL model.

        This function solves for x in:

            four_param_logistic(x, *params) = y_obs

        using a 1-D root finder (Brent's method).

        Parameters
        ----------
        y_obs : float
            Observed OD value.
        params : numpy.ndarray
            Fitted parameters [a, b, c, d].
        x_min : float, optional
            Minimum search bound for the concentration.
        x_max : float, optional
            Maximum search bound for the concentration.

        Returns
        -------
        float
            Estimated concentration corresponding to `y_obs`.
        """
        a, b, c, d = params
        eps = 1e-6

        # Clip the observed OD to be within the modeled range
        y_clipped = min(max(y_obs, d + eps), a - eps)

        def f(x: float) -> float:
            return float(cls.four_param_logistic(x, a, b, c, d) - y_clipped)

        return brentq(f, x_min, x_max)
    
    @classmethod
    def fit_elisa_4pl(
        cls,
        concentrations: np.ndarray,
        od_values: np.ndarray,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Fit the 4PL model to ELISA standard data.

        Parameters
        ----------
        concentrations : numpy.ndarray
            Known concentrations of the standards.
        od_values : numpy.ndarray
            OD readings corresponding to the standards.

        Returns
        -------
        params : numpy.ndarray
            Fitted parameters [a, b, c, d].
        cov : numpy.ndarray
            Parameter covariance matrix from `curve_fit`.
        """
        # Initial parameter guesses:
        a0 = float(od_values.max())
        d0 = float(od_values.min())
        c0 = float(np.median(concentrations))
        b0 = 1.0  # slope

        p0 = [a0, b0, c0, d0]

        params, cov = curve_fit(
            cls.four_param_logistic,
            concentrations,
            od_values,
            p0=p0,
            maxfev=10000,
        )

        return params, cov
    
    @classmethod
    def plot_elisa_curve(
        cls,
        standards: pd.DataFrame,
        samples: pd.DataFrame,
        params: np.ndarray,
        cov: np.ndarray,
        x_min: float = 0.0,
        x_max: float = 250.0,
        output_file: str | None = None,
    ) -> None:
        """
        Plot the fitted ELISA standard curve and sample data using
        a clean, publication-style aesthetic.

        The plot style is inspired by modern scientific figures:
        - Minimal spines (top and right removed).
        - Slightly larger fonts.
        - Orange and teal color palette for the groups.
        - Legend outside or at the edge of the axes.
        """
        # ------------------------------------------------------------------
        # Style configuration
        # ------------------------------------------------------------------
        # Using a similar pallete color
        color_standard = "#F79256"  
        color_sample = "#00B2CA"    

        plt.rcParams.update({
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 11,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 1.2,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.size": 4,
            "ytick.major.size": 4,
            "figure.dpi": 120,
        })

        # Generate smooth x-grid for the fitted curve
        x_plot = np.linspace(x_min, x_max, 200)
        y_fit, y_low, y_up = cls.prediction_band(x_plot, params, cov, alpha=0.05)

        fig, ax = plt.subplots(figsize=(5.5, 3.8))

        # ------------------------------------------------------------------
        # Scatter points
        # ------------------------------------------------------------------
        # Standards: orange filled circles
        ax.scatter(
            standards["Concentration"],
            standards["OD"],
            s=45,
            c=color_standard,
            edgecolors="black",
            linewidths=0.5,
            label="Standard",
            zorder=3,
        )

        # Samples: teal open circles at estimated concentrations
        if "EstimatedConcentration" in samples.columns:
            ax.scatter(
                samples["EstimatedConcentration"],
                samples["OD"],
                s=45,
                facecolors="none",
                edgecolors=color_sample,
                linewidths=1.2,
                label="Measured",
                zorder=3,
            )

        # ------------------------------------------------------------------
        # Fitted curve and confidence band
        # ------------------------------------------------------------------
        ax.plot(x_plot, y_fit, color="black", linewidth=1.5)
        ax.plot(x_plot, y_low, color="black", linestyle="--", linewidth=1.0)
        ax.plot(x_plot, y_up, color="black", linestyle="--", linewidth=1.0)

        # ------------------------------------------------------------------
        # Axes, labels and limits
        # ------------------------------------------------------------------
        ax.set_xlabel("[ng/µL]")
        ax.set_ylabel("OD")
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(0.0, 1.1)

        # Make left and bottom spines a bit thicker
        for spine in ["bottom", "left"]:
            ax.spines[spine].set_linewidth(1.4)

        # Light grid on y-axis (optional, like many modern figures)
        ax.grid(axis="y", alpha=0.2)

        # Legend similar to panel B of your figure
        leg = ax.legend(
            frameon=False,
            loc="upper right",
            borderpad=0.3,
            handlelength=1.5,
        )

        ax.set_title("ELISA (450 nm)")

        fig.tight_layout()

        if output_file is not None:
            fig.savefig(output_file, dpi=300, bbox_inches="tight")

        plt.show()
