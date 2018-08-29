# Classical Optimization

Classical Optimization Algorithms

## Algorithms

### Robust Search

Example:

```julia
# objective function
f(x) = x^2 + 54 / x

a, b, n = 0, 5, 10
robustSearch(f, a, b, n)
```

Output: `(2.5, 3.5)`
