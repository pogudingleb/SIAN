# Example 12 in the paper "Global Identifiability of Differential Models", taken from
# Demignot, S., D., D., 
# Effect of prosthetic sugar groups on the pharmacokinetics of glucose-oxidase
read("IdentifiabilityODE.mpl");

sigma := [
diff(x1(t), t) = 0,
diff(x2(t), t) = 0,
y(t) = x1(t) + a,
y2(t) = x2(t) + b
];

IdentifiabilityODE(sigma, GetParameters(sigma), known_initial_values=[x1(0)]);
