causal-model-constructor
========================

takes files listing speakers and tuples of the form "cause \ relation \ effect" and creates a dictionary of causal models, where the keys are the speakers' names and the values are directed, labeled graphs.

An input file must be of the form:

    Mr. MITCHEL
    Taxes / - / People are working hard
    People are working hard / + / Firms are making profits
    Firms are making profits / + / US GDP
    US GDP / - / Unemployement
    Mr. KENNEDY
    ...
    
The output is a python dictionary that maps from names, such as "Mr. MITCHEL", to directed labeled graphs, where the nodes are concepts, such as "Taxes" and "People are working hard" and the edges are weighted arcs that represent the believed causal relations between the concepts. 


    

