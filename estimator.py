import numpy as np
from scipy.signal import correlate # type: ignore
from pandas.plotting import autocorrelation_plot # type: ignore
from scipy.stats import chi2 # type: ignore
import matlab.engine # type: ignore

class Estimator:
    def __init__(self) -> None:
        self.eng = matlab.engine.start_matlab()

    def calculate_autocovariance_and_spectral_density(self, e, τmax, ω):
        N = len(e)  # Length of the input signal

        # Initialize autocovariance and spectral density arrays
        λe = np.zeros(2 * τmax + 1)
        ϕe = np.zeros(len(ω))
        
        # Calculate autocovariance
        for τ in range(-τmax, τmax + 1):
            if τ >= 0:
                λe[τ + τmax] = np.sum(e[:N - τ] * e[τ:]) / (N - τ)
            else:
                λe[τ + τmax] = np.sum(e[-τ:N] * e[:N + τ]) / N

        # λe = np.correlate(e - np.mean(e), e - np.mean(e), mode='full') # cross-covariance
        
        # Calculate spectral density
        for k in range(len(ω)):
            ϕe[k] = (1/N) * np.abs(np.sum(e * np.exp(-1j * ω[k] * np.arange(N))))
        
        # Create lag vector τ
        τ = np.arange(-τmax, τmax + 1)
        
        return λe, τ, ϕe
    
    def auto_cov_spectra(self, e, max_lag, w):
        lags = np.arange(-max_lag, max_lag + 1)
        lambda_ = correlate(e, e, mode='full') / len(e)  # Autocovariance calculation
        n = len(w)
        phi = np.zeros(n)

        for k in range(n):
            phi[k] = np.sum(np.cos(w[k] * lags) * lambda_)

        return lambda_, lags, phi

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

