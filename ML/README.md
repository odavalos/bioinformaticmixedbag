# Machine Learning

## Interpretable ML

If you want to build interpretable models read this book by Christoph Molnar: [Interpretable Machine Learning](https://christophm.github.io/interpretable-ml-book/)
- Provides advantages and disadvantages to a wide range of interpretability methods


#### Random Forest

Feature Importance
- Do not use feature importance from the random forest model directly!
  - https://explained.ai/rf-importance/index.html

- Calculate the permutation importance
  - [sklearn permutation importance](https://scikit-learn.org/stable/modules/permutation_importance.html)
    - If your dataset is large this becomes very computationally expensive
  - [rfpimp](https://github.com/parrt/random-forest-importances) is a great alternative
- Calculate [Lime (local interpretable model-agnostic explanations)](https://github.com/marcotcr/lime)
- Calculate [Shapley values](https://github.com/slundberg/shap)
  - [Original Paper here](http://papers.nips.cc/paper/7062-a-unified-approach-to-interpreting-model-predictions)


