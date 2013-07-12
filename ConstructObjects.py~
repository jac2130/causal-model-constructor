import string
from GraphWorld import GraphWorld, GraphCanvas, Layout, CircleLayout, RandomLayout
import DirectedGraphWorld, nltk

from DirectedGraph import SignedArc, DirectedGraph

from Graph import Vertex, Graph

class Sgn(object):
    """A Sgn is the sign of causation in a causal graph."""

    def __init__(self, label='0'):
        self.label = label
        if self.label=='0':
            self.color = "blue"
        elif self.label=='+':
            self.color="black"
        elif self.label =='-':
            self.color="red"
        else: self.color="orange"


    def __repr__(self):
        """Returns a string representation of this object that can
        be evaluated as a Python expression."""
        return 'Sgn(%s)' % repr(self.label)
    __str__ = __repr__
    """The str and repr forms of this object are the same."""


class CausalGraph(DirectedGraph):

    def __init__(self,vs=set(),es=set()):
        """
        Creates a new directed graph.
        @Args:
            vs, a list of Vertices
            es, a list of Arcs
        @Returns:
            None
        """
        self.reverse_graph = {} #keeps a map of in-vertices
        for v in vs:
            self.add_vertex(v)
        for e in es:
            self.add_signed_arc(e)

    def add_signed_arc(self, e):
        """
        Creates an arc FROM v TO w (As a value of v in the internal dictionary).
        Adds same arc TO w FROM v (As a value of w in the reverse dictionary).
        @Args:
            An Arc e
        @Returns:
            None
        """
        v, sgn, w = e
        if v == w:
            raise LoopError('An arc cannot exist from a vertex to itself.')

        self[v][w] = e

        self.reverse_graph[w][v] = e

    add_arc = add_edge=add_signed_arc
    """We only want to add arcs, not edges"""

from collections import defaultdict

def alpnum_cycle():
    num=0
    while True:
        for c in string.uppercase:
            if c=='A': num+=1
            yield c + repr(num)


class Discussion(defaultdict):

    def __init__(self, lines=[], people=[], name=''):
        """
        Creates a Representation of a Congressional Discussion.
        @Args:
            lines, a list of lines from output of the causal sentense parser
            each line must currently be either a name of a person or a person's
            statement of one or more causal beliefs, seperated by
            / AND /. If it is a person's name, the dictionary searches for that name,
            which is then either found or not. If it is not found, a new causal graph (or belief system)
            is created and indexed by the person's name in the dictionary of names, mapping to causal beliefs.
            If the name is found, or if it has been entered in the dictionary, it starts updating that person's
            believes with the following lines untul another name is encountered.
            people, a dictionary mapping from names of senators, to their belief systems (causal graphs)
        @Returns:
            None
        """

        self.name=name
        self.relable=defaultdict(str)
        self.people=set(people + self.keys())
        self.new_names=list(set(people)-set(self.keys())) #New People being added.
        self.process_lines(lines)
        self.concepts=set()

        for name in self.new_names:
            self[name]=CausalGraph()


    def process_lines(self, lines):

        for i in range(len(lines)-1):
            if is_name(lines[i]):
                if len(lines[i].split(' '))>2:
                    if lines[i].split(' ')[2].upper()==lines[i].split(' ')[2]:
                        self.name=lines[i].split('.')[2].lower().lstrip().rstrip()
                    else:
                        self.name=lines[i].split('.')[1].lower().lstrip().rstrip()
                else:
                    self.name=lines[i].split('.')[1].lower().lstrip().rstrip()
                name=self.name
                if name not in self.keys():
                    self[name]=CausalGraph()
                self.add_causal_believes(name, lines[i+1])
            else:
                name=self.name
                self.add_causal_believes(name, lines[i])

    def add_causal_believes(self, name, line, alph=alpnum_cycle()):
        #import pdb
        #pdb.set_trace()
        arcs=[sarc_from_string(chunk) for chunk in and_split(line)]
        vs=sorted(set(self[name].keys()))
        vs_labels=[key.label for key in vs]
        vs_short       =[key for key in vs if len(key.label)==2]
        vs_short_labels=set([key.label for key in vs if len(key.label)==2])
        for arc in arcs:
            x, sgn, y = arc

            if self.relable[x.label]=='':


                #vs_labels.append(x.label)


                self.relable[x.label]=alph.next()
                x.label=self.relable[x.label]


            else:
                x.label=self.relable[x.label]



            if self.relable[y.label]=='':

                self.relable[y.label]=alph.next()
                y.label=self.relable[y.label]




            else:
                y.label=self.relable[y.label]

            if x.label in vs_labels:

                x=[v for v in vs_short if v.label==x.label][0]


            if y.label in vs_labels:
                y=[v for v in vs_short if v.label==y.label][0]

            if x.label not in vs_labels:
                vs.append(x)
                self[name].add_vertex(x)
            if y.label not in vs_labels:
                vs.append(y)
                self[name].add_vertex(y)

            arc1= x, sgn, y

            self[name].add_signed_arc(arc1)

def sarc_from_string(line):
    cause=Vertex(line.split("/")[0].lstrip().rstrip()) #stripping out white space on both ends of the string
    sgn=Sgn(line.split("/")[1].lstrip().rstrip())
    effect=Vertex(line.split('/')[2].lstrip().rstrip())
    return SignedArc([cause, sgn, effect])

def is_name(line):
    line=line.split(' ')
    if (line[0]=='Mr.' or line[0]=='Mrs.') and line[1].upper()==line[1]:
        print line
        return True
    else: return False

def and_split(line):
    return line.split("/ AND /")


def main(script, *args):
    path="/home/johannes/Dropbox/CausalAssertions/"
    import os
    os.chdir(path)
    #This is where the output from DataTemple will go and where this output gets converted to a list of signed arcs:

   #try this for a causal graph that is taking as its input all of the actual output from all of our sample sentences and templates. It is visually hard to interpret:
    #es=[sarc_from_string(line) for text in open('outputs.txt', "r") for line in and_split(text)]

    linez=[line.lstrip().rstrip() for text in open('demo_output.txt', "r").read().split("\n") for line in and_split(text)]
#Making sure that concepts that have been considered to be unique are uniquely represented in the causal graph:

    discussion =Discussion(linez)
    import CausalGraphWorld
    key=discussion.keys()[1]
    graph=discussion[key]
    show_causal_graph(graph, write_postscript=True, filename=key + '.ps')
    key_table= sorted([(value, key) for key, value in discussion.relable.items()])
    for value, key in key_table:
        print "%s = %s" % (value, key) #Print out what the new values mean.
    #os.system('ps2pdf -dDEVICEWIDTHPOINTS=600 -dDEVICEHEIGHTPOINTS=600' + ' ' + key +'.ps') #there is some issue that I still have to resolve.
    #os.system("pdfcrop -margins 10 " + key + '.pdf' + ' ' + key + '.pdf')
    #print arc_c
    #print arc_b
    #print vs
    ###print key
    print discussion['fleming'].keys()
    print discussion.keys( )
    #print discussion[''] why is Mrs. MCMORRIS RODGERS. not entered into the discussion object?

if __name__ == '__main__':
    import sys
    main(*sys.argv)
