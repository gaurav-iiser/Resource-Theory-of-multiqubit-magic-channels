# Resource-Theory-of-multiqubit-magic-channels

This repository contains two codes:
1. conv_hull: This code allows the user to easily create stabilizer convex hulls given the Bloch vector correspoding to a particular state. One can add as many convex polytope. This code helps to visualize whether a convex polytope corresponding to one state lies inside the convex polytope of another state. If it does, then we can infer that we can convert the second state to the first state using completely stabilizer preserving operations (CSPOs), else we cannot.
2.  simple_interconversion_qubits: This code takes the Bloch vectors of the initial and target qubits and outputs whether the initial state can or cannot be converted to the target state.
