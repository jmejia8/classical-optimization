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

![Robust Search](https://www.candaana.com/report/wp-content/uploads/2018/08/rs-ejemplo.gif)
