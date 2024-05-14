import matplotlib.pyplot as plt
plt.style.use('science')

x = [0, 1, 2, 3]
y1 = [1.0, 0.6360, 0.5813, 0.5861]
y1_err = [0.0, 0.0001, 0.0001, 0.0016]
y2 = [1.0, 0.6704, 0.5921, 0.6021]
y2_err = [0.0, 0.0001, 0.0002, 0.0006]

y1_err =[e*4.0 for e in y1_err]
y2_err =[e*4.0 for e in y2_err]

colors = ['blue', 'red', 'green', 'black', 'purple', 'brown', 'orange', 'gray']

# Set the figure size and dpi
fig = plt.figure(figsize=(6, 4), dpi=100)

# Plot y1 with error bars
plt.errorbar(x, y1, yerr=y1_err, fmt='o-', color='blue', capsize=5, markersize=5, label='$G_0W_0$')

# Plot y2 with error bars
plt.errorbar(x, y2, yerr=y2_err, fmt='o-', color='red', capsize=5, markersize=5, label='$G_0W_{KO}$')

# Add grid lines
plt.grid(linestyle='--', alpha=0.5)

# Add legend, labels, and title
plt.legend(fontsize=16)
plt.xlabel('order', fontsize=16)
plt.ylabel('$z$', fontsize=16)

# Set the x- and y-axis limits
plt.xlim(-0.5, 3.5)
plt.ylim(0.5, 1.1)

# Set the tick label font size
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)

# Remove the top and right spines
# Save the figure as a PNG file
plt.savefig('convergence.pdf', bbox_inches='tight', dpi=300)

# Display the plot
plt.show()

