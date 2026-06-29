import matplotlib.pyplot as plt

# Data extracted from the computational analysis
compounds = ['LiCl', 'NaCl', 'KCl', 'CsCl']
vd2_values = [361.90, 329.06, 318.69, 313.89]

# Initialize the figure size
plt.figure(figsize=(8, 6))

# Plotting the line graph
# marker='o' creates the circled points
# markersize and markerfacecolor visually highlight the exact data nodes
plt.plot(compounds, vd2_values, marker='o', markersize=10, markerfacecolor='red', 
         linestyle='-', color='blue', linewidth=2.5, label='vd2 Center Frequency')

# Loop to annotate each circled point with its exact numerical value
for i, v in enumerate(vd2_values):
    plt.text(i, v + 3, f"{v:.2f}", ha='center', va='bottom', fontweight='bold', fontsize=10)

# Formatting the graph with labels, title, and grid
plt.title('In-Plane Libration Center across Alkali Chlorides', fontsize=14)
plt.xlabel('1M Solute Compound', fontsize=12)
plt.ylabel('Center Frequency ($cm^{-1}$)', fontsize=12)
plt.grid(True, linestyle='--', alpha=0.6)
plt.ylim(300, 380) # Adjusting y-axis limits to accommodate the text annotations
plt.legend()

# Display the plot in a window instead of saving it
plt.show()