# Simple 2D mass-spring differentiable simulation

```
//////////////////////////
--------------------------
  (0,0) \     //////  
         \    ------
    l1,k1 \     /  (1,0.5)
           \  /  l2,k2
            O   m
         (0.5, 1)  
```

For the given system, the script calculates spring constants that would result in and eqauilibrium positions of the point mass to be `target_position` using the [differentiable simulation concept from the Siggraph course](https://dl.acm.org/doi/10.1145/3476117.3483433).

Feel free to change the parameters in the script and it should still find an optimal solution.
