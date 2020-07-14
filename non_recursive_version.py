import math
import time
import numpy as np
import matplotlib.pyplot as plt
from functions import *
import multiprocessing as mp
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *

hist_all_nodes = []
hist_remaining_nodes = []

hist_time_for_method = []
hist_time_for_tarjan = []

class BaseGraph:
    def __init__(self, graph_dict=None):
        """ initializes a graph object 
            If no dictionary or None is given, 
            an empty dictionary will be used
        """
        self.Time = 0
        self.scc_field = {}

        if graph_dict == None:
            graph_dict = {}
        self.__graph_dict = graph_dict

    def vertices(self):
        """ returns the vertices of a graph """
        return list(self.__graph_dict.keys())

    def edges(self):
        """ returns the edges of a graph """
        return self.__generate_edges()

    def add_vertex(self, vertex):
        """ If the vertex "vertex" is not in 
            self.__graph_dict, a key "vertex" with an empty
            list as a value is added to the dictionary. 
            Otherwise nothing has to be done. 
        """
        if vertex not in self.__graph_dict:
            self.__graph_dict[vertex] = []

    def has_vertex(self, vertex):
        if vertex in self.__graph_dict:
            return 1
        return 0

    def number_of_vertices(self):
        return len(self.vertices())

    def has_edge(self, edge):
        key = edge[0]
        value = edge[1]
        return value in self.__graph_dict[key]

    def add_edge(self, edge):
        """ assumes that edge is of type set, tuple or list; 
            between two vertices can be multiple edges! 
        """
        if(self.has_edge(edge) == 0):
            self.__graph_dict[edge[0]].append(edge[1])

    def number_of_edges(self):
        number_of_edges = 0
        for key in self.__graph_dict:
            number_of_edges += len(self.__graph_dict[key])
        return number_of_edges

    def __generate_edges(self):
        """ A static method generating the edges of the 
            graph "graph". Edges are represented as sets 
            with one (a loop back to the vertex) or two 
            vertices 
        """
        edges = []
        for vertex in self.__graph_dict:
            for neighbour in self.__graph_dict[vertex]:
                edges.append({vertex, neighbour})
        return edges

    def SCC(self):
        next_id = 0 # next_id index.
        length = max(self.vertices()) + 1
        index = [None] * length
        lowlink = [None] * length
        onstack = [False] * length
        stack = []
        next_idgroup = 0 # next_id SCC ID.
        groups = [] # SCCs: list of vertices.
        groupid = {} # Map from vertex to SCC ID.
        for v in self.vertices():
            if index[v] == None:
                self.sconnect(v, next_id, next_idgroup, index, lowlink, onstack, stack, groups, groupid)

    def sconnect(self, v, next_id, next_idgroup, index, lowlink, onstack, stack, groups, groupid):
        work = [(v, 0)] # NEW: Recursion stack.
        #k = 0
        while work:
            v, i = work[-1] # i is next_id successor to process.
            del work[-1]
            if i == 0: # When first visiting a vertex:
                index[v] = next_id
                lowlink[v] = next_id
                next_id += 1
                stack.append(v)
                onstack[v] = True
            recurse = False
            for j in range(len(self.__graph_dict[v])):
                w = self.__graph_dict[v][j]
                if index[w] == None:
                    # CHANGED: Add w to recursion stack.
                    work.append((v, j+1))
                    work.append((w, 0))
                    recurse = True
                    break
                elif onstack[w]:
                    lowlink[v] = min(lowlink[v], index[w])
            if recurse: continue # NEW
            k = 0
            if index[v] == lowlink[v]:
                com = []
                while True:
                    w = stack[-1]
                    del stack[-1]
                    onstack[w] = False
                    com.append(w)
                    if(w != v or k  or self.has_edge([w, w])):
                        self.scc_field[w] = 1
                        k = 1
                    groupid[w] = next_idgroup
                    if w == v: break
                groups.append(com)
                next_idgroup += 1
            if work: # NEW: v was recursively visited.
                w = v
                v, _ = work[-1]
                lowlink[v] = min(lowlink[v], lowlink[w])

    def delete_scc_components(self):
        start = timerStart()
        self.SCC()
        print("------------ Time for Tarjan: ", timeFinish(start))
        # print(self.__graph_dict)
        start = timerStart()
        for vertex in self.vertices():
            if(vertex not in self.scc_field):
                self.__graph_dict.pop(vertex)
        # print("------------", self.__graph_dict)
        print("------------ Time for deletion: ", timeFinish(start))

    def __str__(self):
        res = str(self.__graph_dict)
        res += "\nvertices: "
        for k in self.__graph_dict:
            res += str(k) + " "
        return res

class FinalGraph():

    def __init__(self, givenBoundaries, cell_size, cell_count=0, lengthOfSide=[0, 0], Time=0):
        self.graph = BaseGraph()
        self.boundaries = givenBoundaries
        self.cell_size = cell_size
        self.cell_count = cell_count
        self.horizontal_length = lengthOfSide[0]
        self.vertical_length = lengthOfSide[1]

    def calculate_cell_count(self):
        self.cell_count = abs(((self.boundaries[1][0] - self.boundaries[0][0]) / self.cell_size) * ((self.boundaries[1][1] - self.boundaries[0][1]) / self.cell_size))

    def calculate_sides_of_area(self):
        self.horizontal_length = int(abs((self.boundaries[1][0] - self.boundaries[0][0]) / self.cell_size))
        self.vertical_length = int(abs((self.boundaries[1][1] - self.boundaries[0][1]) / self.cell_size))

    def create_cells(self):
        for row in range(self.vertical_length):
            for col in range(self.horizontal_length):
                self.graph.add_vertex(col + row * self.horizontal_length)

    def get_local_coordinates_of_cell_from_number(self, number_of_cell):
        return [int(math.floor(number_of_cell / self.horizontal_length)), number_of_cell % self.vertical_length]
    
    def get_local_coordinates_of_cell_from_global(self, coord):
        # col = math.floor(
        #     (coord[0] * self.horizontal_length / (self.boundaries[1][0] - self.boundaries[0][0])) + self.horizontal_length / 2)
        # row = math.floor(
        #     (coord[1] * self.vertical_length / (self.boundaries[0][1] - self.boundaries[1][1])) + self.vertical_length / 2)
        # row = math.fabs(row - self.vertical_length + 1)

        col = math.floor((coord[0] - self.boundaries[0][0]) / self.cell_size)
        row = math.floor((self.boundaries[0][1]- coord[1]) / self.cell_size)

        return [row, col]

    def get_local_coordinates_of_cell_from_global_advanced(self, coord):
        all_local_coord = []
        col = math.floor((coord[0] - self.boundaries[0][0]) / self.cell_size)
        row = math.floor((self.boundaries[0][1]- coord[1]) / self.cell_size)
        all_local_coord.append([row, col])
        if((coord[0] - self.boundaries[0][0]) / self.cell_size == col):
            new_col = col - 1
            if(new_col >= 0):
                all_local_coord.append([row, new_col])
        if((self.boundaries[0][1]- coord[1]) / self.cell_size == row):
            new_row = row - 1
            if(new_row >= 0):
                all_local_coord.append([new_row, col])
        if((coord[0] - self.boundaries[0][0]) / self.cell_size == col and (self.boundaries[0][1]- coord[1]) / self.cell_size == row):
            new_col = col - 1
            new_row = row - 1
            if(new_col >= 0 and new_row >=0):
                all_local_coord.append([new_row, new_col])
        return all_local_coord

    
    def get_global_coordinates_of_cell_from_number(self, number_of_cell):
        x = self.boundaries[0][0] + self.cell_size * self.get_local_coordinates_of_cell_from_number(number_of_cell)[1] 
        y = self.boundaries[1][0] - self.cell_size * self.get_local_coordinates_of_cell_from_number(number_of_cell)[0]
        return [x, y]
    
    def get_number_of_cell_from_local(self, coord):
        return coord[1] + self.horizontal_length * coord[0]
    
    def get_number_of_cell_from_global(self, coord):
        return int(self.get_number_of_cell_from_local(self.get_local_coordinates_of_cell_from_global(coord)))

    def get_number_of_cell_from_global_advanced(self, coord):
        all_numbers = []
        all_local_coord = self.get_local_coordinates_of_cell_from_global_advanced(coord)
        for local_coord in all_local_coord:
            all_numbers.append(self.get_number_of_cell_from_local(local_coord))
        return all_numbers
    
    def check_if_point_inside_boundaries(self, coord):
        if (coord[0] < self.boundaries[0][0]):
            return 0
        elif (coord[0] > self.boundaries[1][0] + self.cell_size):
            return 0
        if (coord[1] > self.boundaries[0][1]):
            return 0
        elif (coord[1] < self.boundaries[1][1] - self.cell_size):
            return 0
        return 1

    def clip_point_to_boundaries(self, coord):
        x = coord[0]
        y = coord[1]
        if (x == self.boundaries[0][0]):
            x = self.boundaries[0][0]
        elif (x == self.boundaries[1][0]):
            x = self.boundaries[1][0] - self.cell_size

        if (y == self.boundaries[0][1]):
            y = self.boundaries[0][1]
        elif (y == self.boundaries[1][1]):
            y = self.boundaries[1][1] + self.cell_size
        return [x, y]

    def fill_nodes(self):
        self.calculate_cell_count()
        self.calculate_sides_of_area()
        self.create_cells()

    def hit_cell(self, coord):
        coord = self.clip_point_to_boundaries(coord)
        if(self.check_if_point_inside_boundaries(coord)):
            number_of_hit_cell = self.get_number_of_cell_from_global(coord)
            if self.graph.has_vertex(number_of_hit_cell):
                return number_of_hit_cell
        return -1

    def hit_cell_advanced(self, coord):
        hit_cells = []
        coord = self.clip_point_to_boundaries(coord)
        if(self.check_if_point_inside_boundaries(coord)):
            numbers_of_hit_cell = self.get_number_of_cell_from_global_advanced(coord)
            for number in numbers_of_hit_cell:
                if self.graph.has_vertex(number):
                    hit_cells.append(number)
        if(len(hit_cells) != 0):
            return hit_cells
        else:
            return -1

    def point_method(self, increment = 1):
        for vertex in self.graph.vertices():
            for j in range(increment):
                y = self.get_global_coordinates_of_cell_from_number(vertex)[1] - float(j * self.cell_size) / increment
                for i in range(increment):
                    x = self.get_global_coordinates_of_cell_from_number(vertex)[0] + float(i * self.cell_size) / increment              
                    number_of_hit_cell = self.hit_cell(f([x, y]))
                    if number_of_hit_cell != - 1:
                        self.graph.add_edge([vertex, int(number_of_hit_cell)])
    def point_method_advanced(self, increment = 1):
        for vertex in self.graph.vertices():
            for j in range(increment):
                y = self.get_global_coordinates_of_cell_from_number(vertex)[1] - float(j * self.cell_size) / increment
                for i in range(increment):
                    x = self.get_global_coordinates_of_cell_from_number(vertex)[0] + float(i * self.cell_size) / increment              
                    numbers_of_hit_cell = self.hit_cell_advanced(f([x, y]))
                    if numbers_of_hit_cell != - 1:
                        for number in numbers_of_hit_cell:
                            self.graph.add_edge([vertex, int(number)])


    def hit_area(self, vertex, top_left_number, bot_right_number, top_right_number, bot_left_number):
        all_local_coord_row = []
        all_local_coord_col = []
        
        all_local_coord_row.append(self.get_local_coordinates_of_cell_from_number(top_left_number)[0])
        all_local_coord_row.append(self.get_local_coordinates_of_cell_from_number(top_right_number)[0])
        all_local_coord_row.append(self.get_local_coordinates_of_cell_from_number(bot_left_number)[0])
        all_local_coord_row.append(self.get_local_coordinates_of_cell_from_number(bot_right_number)[0])

        all_local_coord_col.append(self.get_local_coordinates_of_cell_from_number(top_left_number)[1])
        all_local_coord_col.append(self.get_local_coordinates_of_cell_from_number(top_right_number)[1])
        all_local_coord_col.append(self.get_local_coordinates_of_cell_from_number(bot_left_number)[1])
        all_local_coord_col.append(self.get_local_coordinates_of_cell_from_number(bot_right_number)[1])

        col_min = min(all_local_coord_col)
        col_max = max(all_local_coord_col)
        row_min = min(all_local_coord_row)
        row_max = max(all_local_coord_row)

        low_bound = 1
        high_bound = 2

        for row in range(int(row_min) - low_bound, int(row_max) + high_bound):
            for col in range(int(col_min) - low_bound, int(col_max) + high_bound):
                number_of_cell = self.get_number_of_cell_from_local([row, col])
                if self.graph.has_vertex(number_of_cell):
                    self.graph.add_edge([vertex, number_of_cell])

    def linear_method(self):
        top_left_corner_coord = [0, 0]
        bottom_right_corner_coord = [0, 0]
        for vertex in self.graph.vertices():
            top_left_corner_coord = self.get_global_coordinates_of_cell_from_number(vertex)
            bottom_left_corner_coord = [self.get_global_coordinates_of_cell_from_number(vertex)[0], 
                                         self.get_global_coordinates_of_cell_from_number(vertex)[1] - self.cell_size]
            top_right_corner_coord = [self.get_global_coordinates_of_cell_from_number(vertex)[0] + self.cell_size, 
                                      self.get_global_coordinates_of_cell_from_number(vertex)[1]]
            bottom_right_corner_coord = [self.get_global_coordinates_of_cell_from_number(vertex)[0] + self.cell_size, 
                                         self.get_global_coordinates_of_cell_from_number(vertex)[1] - self.cell_size] 

            hit_coord_top_left = f(top_left_corner_coord)
            hit_coord_top_right = f(top_right_corner_coord)
            hit_coord_bot_left = f(bottom_left_corner_coord)
            hit_coord_bot_right = f(bottom_right_corner_coord)

            #EXPERIMENTAL
            amount_of_mapping = 0
            for i in range(amount_of_mapping):
                hit_coord_top_left = f(hit_coord_top_left)
                hit_coord_bot_right = f(hit_coord_bot_right)
                hit_coord_top_right = f(hit_coord_top_right)
                hit_coord_bot_left = f(hit_coord_bot_left)

            number_of_hit_cell_top_left = self.hit_cell(hit_coord_top_left)
            number_of_hit_cell_top_right = self.hit_cell(hit_coord_top_right)
            number_of_hit_cell_bot_left = self.hit_cell(hit_coord_bot_left)
            number_of_hit_cell_bot_right = self.hit_cell(hit_coord_bot_right)

            #DEBUG
            
            # if(hit_coord_top_left[0] > 1.1 and hit_coord_top_left[0] < 1.7) and (hit_coord_top_left[1] > 0.9 and hit_coord_top_left[1] < 1.4):
            #     print("Vertex that maps to that point: ", vertex)
            #     print("Coordinates of that vertex: ", top_left_corner_coord)
            #     print("Vertex that is that point: ", number_of_hit_cell_top_left)
            #     print("Hit coordinates: ", hit_coord_top_left)

            # if(hit_coord_bot_right[0] > 1.1 and hit_coord_bot_right[0] < 1.7) and (hit_coord_bot_right[1] > 0.9 and hit_coord_bot_right[1] < 1.4):
            #     print("Vertex that maps to that point: ", vertex)
            #     print("Coordinates of that vertex: ", bottom_right_corner_coord)
            #     print("Vertex that is that point: ", number_of_hit_cell_bot_right)
            #     print("Hit coordinates: ", hit_coord_bot_right)
            # print("Vertex to map: ", vertex)
            # print("Coord to map: ", top_left_corner_coord, bottom_right_corner_coord, top_right_corner_coord, bottom_left_corner_coord)
            # print("Mapped vertices: ", number_of_hit_cell_top_left, number_of_hit_cell_bot_right, number_of_hit_cell_top_right, number_of_hit_cell_bot_left)
            # print("Mapped coord: ", hit_coord_top_left, hit_coord_bot_right, hit_coord_top_right, hit_coord_bot_left)

            
            if(number_of_hit_cell_bot_right != -1 and number_of_hit_cell_top_left != -1 and number_of_hit_cell_bot_left != -1 and number_of_hit_cell_top_right != -1):
                self.hit_area(vertex, number_of_hit_cell_top_left, number_of_hit_cell_bot_right, number_of_hit_cell_top_right, number_of_hit_cell_bot_left)
            else:
                if(number_of_hit_cell_bot_right != -1):
                    if(self.graph.has_vertex(number_of_hit_cell_bot_right)):
                        self.graph.add_edge([vertex, number_of_hit_cell_bot_right])
                if(number_of_hit_cell_top_left != -1):
                    if(self.graph.has_vertex(number_of_hit_cell_top_left)):
                        self.graph.add_edge([vertex, number_of_hit_cell_top_left])
                if(number_of_hit_cell_bot_left != -1):
                    if(self.graph.has_vertex(number_of_hit_cell_bot_left)):
                        self.graph.add_edge([vertex, number_of_hit_cell_bot_left])
                if(number_of_hit_cell_top_right != -1):
                    if(self.graph.has_vertex(number_of_hit_cell_top_right)):
                        self.graph.add_edge([vertex, number_of_hit_cell_top_right])


    def prepare_new_graph(self):
        new_graph = FinalGraph(self.boundaries, self.cell_size / 2)
        new_graph.cell_count = self.cell_count * 4
        new_graph.horizontal_length = self.horizontal_length * 2
        new_graph.vertical_length = self.vertical_length * 2
        return new_graph

    def get_four_new_coordinates(self, coord):
        eps = self.cell_size / 16
        return [[coord[0] + eps, coord[1] - eps], 
               [coord[0] + self.cell_size + eps, coord[1] - eps], 
               [coord[0] + eps, coord[1] - self.cell_size - eps], 
               [coord[0] + self.cell_size + eps, coord[1] - self.cell_size - eps]]

    def divide_scc_cells(self):
        print("Number of vertices: ", len(self.graph.vertices()))
        graph_after_division = self.prepare_new_graph()
        for vertex in self.graph.vertices():
            new_cell_coordinates = graph_after_division.get_four_new_coordinates(self.get_global_coordinates_of_cell_from_number(vertex))
            for coord in new_cell_coordinates:
                graph_after_division.graph.add_vertex(int(graph_after_division.get_number_of_cell_from_global(coord)))
        return graph_after_division

    def full_alghorithm(self):
        print("Current step: ", amount_of_steps_taken)
        time = timerStart()
        if(self.cell_size == cell_size):
            self.fill_nodes()
            print("Time to fill nodes: ", timeFinish(time))
        time = timerStart()
        amount_of_points = 4
        #self.point_method(amount_of_points)
        #self.point_method_advanced(amount_of_points)
        self.linear_method()
        print("Time for linear method: ", timeFinish(time))
        hist_time_for_method.append(timeFinish(time))
        time = timerStart()
        hist_all_nodes.append(self.graph.number_of_vertices())
        self.graph.delete_scc_components()
        hist_time_for_tarjan.append(timeFinish(time))
        print("Time to delete scc components: ", timeFinish(time))

        #DEBUG
        # k = 0
        # print("Vertices with edges to itself")
        # for v in self.graph.vertices():
        #     if self.graph.has_edge([v, v]):
        #         print(v)
        #         k += 1
        # print("Amount: ", k)
        # print("Relative Amount: ", k / len(self.graph.vertices()))
        hist_remaining_nodes.append(self.graph.number_of_vertices())
        print("Number of vertices: ", self.graph.number_of_vertices())
        print("Number of edges: ", self.graph.number_of_edges())

        if amount_of_steps_taken < amount_of_steps - 1:
            time = timerStart()
            final_graph = self.divide_scc_cells()
            print("Time to divide cells: ", timeFinish(time))
            return final_graph
        else:
            return self


    def __str__(self):
        res = "vertices: "
        res += str(self.graph.vertices())
        res += "\nnumber of vertices: "
        res += str(len(self.graph.vertices()))
        res += "\nedges: "
        res += str(self.graph.edges())
        res += "\nnumber of cells: " + str(self.cell_count)
        res += "\nlengths of sides: " + str(self.horizontal_length) + ", " + str(self.vertical_length)
        return res

def f_with_bounds(coord):
    return 

def timerStart():
    return time.time()
def timeFinish(mTime):
    return time.time()-mTime

def drawLine():
    glBegin(GL_LINES)
    # glColor3f(1, 0, 0)
    glVertex2f(cBound[0][0], 0)
    glVertex2f(abs(cBound[1][0]), 0)
    glVertex2f(0, cBound[0][1])
    glVertex2f(0, cBound[1][1])
    glEnd()


def drawGraph(G):
    for node in G.graph.vertices():
        glRectd(G.get_global_coordinates_of_cell_from_number(node)[0], G.get_global_coordinates_of_cell_from_number(node)[1],
                G.get_global_coordinates_of_cell_from_number(node)[0] + G.cell_size, G.get_global_coordinates_of_cell_from_number(node)[1] - G.cell_size)


def init():
    glClearColor(1, 1, 1, 1)
    glMatrixMode(GL_PROJECTION)
    glLoadIdentity()

    glOrtho(cBound[0][0], math.fabs(cBound[1][0]), cBound[1][1], math.fabs(cBound[0][1]), 0, 1)
    glColor3f(0, 0, 0)
    glPointSize(3)

def display():
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glColor3f(0, 0, 0)
    glPointSize(5)
    glLineWidth(1)

    drawGraph(graph)

    glColor3f(1, 1, 0)
    glPointSize(5)
    glLineWidth(1)

    drawLine()

    glutSwapBuffers()

def my_mouse_routine(button, state, x, y):
    if(button == GLUT_LEFT_BUTTON):
        if(state == GLUT_UP):
            viewport     = glGetIntegerv(GL_VIEWPORT)
            matrixModelView  = glGetDoublev(GL_MODELVIEW_MATRIX)
            matrixProjection = glGetDoublev(GL_PROJECTION_MATRIX)


            realy = viewport[3] - y - 1
            print ('Coordinates at cursor are', x, viewport[3] - y - 1)
            print ('World coords at z=0 are', gluUnProject(x, realy, 0, matrixModelView, matrixProjection, viewport))
            print ('Number of cell: ', graph.get_number_of_cell_from_global([gluUnProject(x, realy, 0, matrixModelView, matrixProjection, viewport)[0], 
                                                                             gluUnProject(x, realy, 0, matrixModelView, matrixProjection, viewport)[1]]))

cBound = cBoundVDP
cell_size = 1.0
amount_of_steps = 11
# TEST
fulltime = timerStart()
graph = FinalGraph(cBound, cell_size)
for amount_of_steps_taken in range(amount_of_steps):
    graph = graph.full_alghorithm()
print("compl")
print("Full time of completion: ", timeFinish(fulltime))

labels = []
for i in range(len(hist_remaining_nodes)):
    labels.append('Step ' + str(i+1))
x = np.arange(len(hist_remaining_nodes))


fig, ax = plt.subplots()

bar1 = ax.bar(x, hist_all_nodes, color='blue', width=0.25, label='Количество всех вершин')
bar2 = ax.bar(x + 0.25, hist_remaining_nodes, color='pink', width=0.25, label='Количество вершин после удаления')

ax.set_xticks(x)
ax.set_xticklabels(labels)



ax.set_ylabel('Количество вершин')
ax.set_xlabel('Количество всех и оставшихся после удаления вершин на примере системы Хенона')
ax.legend()

def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = rect.get_height()
        ax.annotate('{}'.format(height),
                    xy=(rect.get_x() + rect.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

autolabel(bar1)
autolabel(bar2)
fig.tight_layout()

plt.show()


glutInit()
glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB)
glutInitWindowSize(700, 700)
glutInitWindowPosition(200, 250)
glutCreateWindow("Graph")

init()

glutDisplayFunc(display)

glutMouseFunc( my_mouse_routine );
glutMainLoop()
