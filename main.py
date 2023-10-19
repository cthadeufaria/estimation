import numpy as np
import matplotlib.pyplot as plt
from linear_estimator import LinearEstimator

def run_estimator_1():
    sigma = 4
    i = 2000
    e1 = np.random.normal(0, np.sqrt(sigma), i) # Generate a random signal for testing
    tmax = 20  # Maximum lag
    n = 200
    w1c = np.linspace(-np.pi, np.pi, n)  # Frequency vector
    lambda1c, lag1c, phi1c = est.calculate_autocovariance_and_spectral_density(e1, tmax, w1c)

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

def main():
    run_estimator_1()

if __name__ == '__main__':
    est = LinearEstimator()
    main()