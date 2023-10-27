import numpy as np
import matplotlib.pyplot as plt
from estimator import LinearEstimator


def main():
    sigma = 4
    i = 2000
    e1 = np.random.normal(0, np.sqrt(sigma), i)
    tmax = 20  # Maximum lag
    n = 200
    w1c = np.linspace(-np.pi, np.pi, n)  # Frequency vector
    lambda1c, lag1c, phi1c = est.auto_cov_spectra(
        e1, tmax, w1c, method='MATLAB'
    )

    # Plot autocovariance and spectral density
    plt.figure(figsize=(10, 6))
    plt.subplot(2, 1, 1)
    plt.stem(lag1c, lambda1c)
    plt.xlabel('Lag (τ)')
    plt.ylabel('Autocovariance (λe)')
    plt.title('Autocovariance Sequence')

    plt.subplot(2, 1, 2)
    plt.plot(w1c, phi1c)
    plt.xlabel('Frequency (ω)')
    plt.ylabel('Spectral Density (ϕe)')
    plt.title('Spectral Density')

    plt.tight_layout()
    plt.show()

    # QUESTION: What are the differences in calculating estimated and the theorectical power density spectrum?

    alpha = 0.05
    s, t, k_alpha = est.whiteness_test(e1, tmax, alpha, method='Python')


if __name__ == '__main__':
    est = LinearEstimator()
    main()
