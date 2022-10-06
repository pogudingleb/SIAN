# Example 12 in the paper "Global Identifiability of Differential Models", taken from
# Demignot, S., D., D., 
# Effect of prosthetic sugar groups on the pharmacokinetics of glucose-oxidase
read("IdentifiabilityODE.mpl");

sigma := [
diff(x1(t), t) = a,
diff(x2(t), t) = b,
y(t) = x1(t) + x2(t)
];

IdentifiabilityODE(sigma, GetParameters(sigma), known_initial_values=[x1(0)]);
