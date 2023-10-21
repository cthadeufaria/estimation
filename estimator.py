import numpy as np
from scipy.stats import chi2
import matlab.engine


class Estimator:
    def __init__(self) -> None:
        self.eng = matlab.engine.start_matlab()

    def auto_cov_spectra(self, e, τmax, ω, method=None):
        N = len(e)

        λe = np.zeros(2 * τmax + 1)
        ϕe = np.zeros(len(ω))

        if method == 'MATLAB':
            λe = np.array(
                self.eng.xcov(e, τmax, 'unbiased')[0]
            ).T

            ϕe = self.eng.fft(λe, len(ω))

        else:
            for τ in range(-τmax, τmax + 1):
                if τ >= 0:
                    λe[τ + τmax] = np.sum(e[:N - τ] * e[τ:])  # / (N - τ)
                else:
                    λe[τ + τmax] = np.sum(e[-τ:N] * e[:N + τ])  # / N

            for w in range(len(ω)):
                ϕe[w] = np.abs(
                    np.sum(
                        λe * np.exp(-1j * ω[w] * np.arange(len(λe)))
                    )
                )  # * (1/N)

        τ = np.arange(-τmax, τmax + 1)

        return λe, τ, ϕe

    def whiteness_test(e, τmax, alpha):
        N = len(e)  # Length of the input signal

        # Calculate the normalized autocovariance sequence s(τ)
        s = np.correlate(e, e, mode='full') / (N * np.var(e))

        # Create a lag vector τ
        τ = np.arange(-τmax, τmax + 1)

        # Calculate the Ljung-Box test statistic Q
        Q = N * (N + 2) * np.sum(s**2 / (N - τ))

        # Calculate the critical value (confidence interval limit) kα
        df = τmax
        k_alpha = chi2.ppf(1 - alpha, df)

        # Perform the whiteness test
        is_white = Q <= k_alpha

        return s, τ, k_alpha, is_white
