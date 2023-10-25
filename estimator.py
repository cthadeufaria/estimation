import numpy as np
from scipy.stats import chi2
import matlab.engine
import matplotlib.pyplot as plt


class LinearEstimator:
    def __init__(self) -> None:
        self.eng = matlab.engine.start_matlab()


    def auto_cov_spectra(self, e, tmax, ω, method=None):
        N = len(e)

        λe = np.zeros(2 * tmax + 1)
        ϕe = np.zeros(len(ω))

        if method == 'MATLAB':
            λe = np.array(
                self.eng.xcov(e, tmax, 'unbiased')[0]
            ).T

            ϕe = self.eng.fft(λe, len(ω))

        else:
            for t in range(-tmax, tmax + 1):
                if t >= 0:
                    λe[t + tmax] = np.sum(e[:N - t] * e[t:])  # / (N - t)
                else:
                    λe[t + tmax] = np.sum(e[-t:N] * e[:N + t])  # / N

            for w in range(len(ω)):
                ϕe[w] = np.abs(
                    np.sum(
                        λe * np.exp(-1j * ω[w] * np.arange(len(λe)))
                    )
                )  # * (1/N)

        t = np.arange(-tmax, tmax + 1)

        return λe, t, ϕe


    def whiteness_test_b(self, e, tmax, alpha):
        N = len(e)  # Length of the input signal

        # Calculate the normalized autocovariance sequence s(τ)
        s = np.correlate(e, e, mode='full') / (N * np.var(e))
        # s = np.array(
        #         self.eng.xcov(e, tmax, 'normalized')[0]
        #     ).T

        # Create a lag vector τ
        t = np.arange(-tmax, tmax + 1)

        # Calculate the Ljung-Box test statistic Q
        Q = N * (N + 2) * np.sum(s**2 / (N - t))

        # Calculate the critical value (confidence interval limit) kα
        df = tmax
        k_alpha = chi2.ppf(1 - alpha, df)

        # Perform the whiteness test
        is_white = Q <= k_alpha

        return s, t, k_alpha, is_white


    def whiteness_test(self, e, tmax, alpha, method='MATLAB'):
        N = len(e)  # Length of the input signal

        if method == 'MATLAB':
            s = self.eng.xcov(e, tmax, 'normalized')  # Why the output has length (tmax * 2) + 1?
            s = np.array(s).T

            t = np.arange(-tmax, tmax + 1)

        elif method == 'Python':
            s = np.zeros(tmax + 1)
            lambda_e = np.zeros(tmax + 1)

            for t in range(0, tmax + 1):
                lambda_e[t] = np.sum(
                    e[tmax + 1:N - t] * e[(t + tmax + 1):]
                ) / (N - tmax)
            
            s = lambda_e / lambda_e[0]  # Is xcov's output normalized autocovariance sequence already the aforementioned quocient?

            t = np.arange(0, tmax + 1)

        k_alpha = self.eng.norminv(1 - (alpha / 2))  # What exactly does this functions calculate?
        
        threshold = k_alpha / np.sqrt(N - tmax)

        self.plot_normalized_autocov_seq(s, t, threshold)

        return s, t, k_alpha


    def plot_normalized_autocov_seq(self, s, t, threshold) -> None:
        plt.figure(figsize=(10, 6))
        plt.plot(t, s)
        plt.plot(t, np.ones((len(t),)) * threshold)
        plt.plot(t, np.ones((len(t),)) * -threshold)
        plt.xlabel('Lag (τ)')
        plt.ylabel('Distribution s(τ)')
        plt.title('Distribution s(τ) = λ̂ee(τ) / λ̂ee(0)')

        plt.tight_layout()
        plt.show()
