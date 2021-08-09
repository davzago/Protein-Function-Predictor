- My work here revolves around the Cafa challenge (3CFU project + thesis)
- brief introduction of CAFA challenge 
- what actually is the task of the challenge
  - What is multilabel prediction (multiclass & multioutput)
  - The labels have a connection: the ontology (wikipida definition?)
  - the ontology is a DAG (example with subontology) 

**Evaluator**

- discuss my implementation
  - what are the inputs
    - discuss the fact that the program can take one or multiple predictions  
    - discuss the fact that the implementation automatically splits the ontologies based on namespace
  - representation of dag (adjacency list or matrix)
  - metrics
  - what are the outputs and show results 

**Thesis**

- My project is to steal from competitor
- discuss components
  - Naive
  - Blast
  - Logistic regression components
    - what is logistic regression
    - omission of LR proFET and K-mer
    - How InterPRO works 
      - transforming a multilabel problem in a classification problem
  - show key component (show some titles of tree boosting winning competitions)
    - what is a tree
    - what is boosting
    - how LambdaMART works
- give some present results 

**END**
$$
F = \frac{2 \cdot pr \cdot rc}{pr + rc}
$$

$$
S = \sqrt{ru^2\cdot mi^2}
$$

$$
S_{Naive}(G_i,P_j) = \frac{N_{G_i}}{N_D}
$$

$$
S_{Blast}(G_i,P_j) = \frac{\sum_{p \in H_j}I(G_i,p)\cdot B(P_j,p)}{\sum_{p \in H_j}B(P_j,p)}
$$

$$
P(Y=1|X) = \beta_0+\beta_1X
$$

$$
P(Y=1|X) = \frac{e^{\beta_0 + \beta_1X}}{1+e^{\beta_0 + \beta_1X}}
$$

