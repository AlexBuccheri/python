from sklearn.datasets import fetch_openml

fetch = False
if fetch:
    mnist = fetch_openml('mnist_784', version=1)

# Note, scikit learn will cache these data sets in $HOME/scikit_learn_data
print(mnist.keys())