import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
# Corrected data based on user feedback
# Regressao L11/Le do livro do Pope 2000, gráfico 6.24

def power_law(x, a, b, c):
    return a * (x**b) + c

corrected_data_points = {
    "Re_lambda": [30, 70, 300, 1000, 3000, 10000],  # log scale (10^1 to 10^4)
    "L11_L":     [0.8, 0.6, 0.5, 0.43, 0.43, 0.43]
}

# Convert corrected data to a pandas DataFrame
corrected_data = pd.DataFrame(corrected_data_points)

# Redo the curve fitting with corrected data
x_data_corrected = np.log10(corrected_data["Re_lambda"])
y_data_corrected = corrected_data["L11_L"]

# Fit the model to the corrected data
params_corrected, _ = curve_fit(power_law, corrected_data["Re_lambda"], y_data_corrected, p0=[1, -0.1, 0.3])

# Generate smooth curve for plotting with corrected data
x_fit_corrected = np.logspace(1, 4, 500)  # Range of Re_lambda (10^1 to 10^4)
y_fit_corrected = power_law(x_fit_corrected, *params_corrected)

# Plot the corrected data and the fitted model
plt.figure(figsize=(8, 6))
plt.scatter(corrected_data["Re_lambda"], y_data_corrected, label="Corrected data points", color="red")
plt.plot(x_fit_corrected, y_fit_corrected, label="Fitted model (corrected)", color="blue")
plt.xscale("log")
plt.xlabel(r"$R_\lambda$ (log scale)", fontsize=12)
plt.ylabel(r"$L_{11}/L$", fontsize=12)
plt.title("Model Fitting (Corrected): $L_{11}/L$ vs $R_\lambda$", fontsize=14)
plt.legend()
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.show()

# Display corrected model parameters
params_corrected
