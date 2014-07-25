"""

PyVisML: A complete python API for interaction with VisANT VisML xml files


Copyright 2013 Michael Seiler
Boston University
miseiler@gmail.com

This file is part of PyVisML.

PyVisML is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PyVisML is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PyVisML.  If not, see <http://www.gnu.org/licenses/>.


"""

import xml.etree.ElementTree as xml
import numpy as N
from warnings import warn

import xml.dom.minidom as minidom
import StringIO, os, sys

VERSION = '1.35'
DEFAULT_SPECIES = 'uno'

def determine_path():
    """Borrowed from wxglade.py"""
    try:
        root = __file__
        if os.path.islink(root):
            root = os.path.realpath(root)
        return os.path.dirname(os.path.abspath(root))
    except:
        print "I'm sorry, but something is wrong."
        print "There is no __file__ variable. Please contact the author."
        sys.exit()

DEFAULT_METHODS_FILE = os.path.join(determine_path(),'data','default_methods.xml')

try:
    from matplotlib.colors import ColorConverter as cc
    COLORCONVERTER = cc()
except:
    COLORCONVERTER = None


def create_empty_visml_tree(default_methods_file=DEFAULT_METHODS_FILE, **attrib):
    r = VisMLTree()

    if 'ver' not in attrib:
        attrib.update({'ver': VERSION})

    if 'species' not in attrib:
        attrib.update({'species': DEFAULT_SPECIES})

    if 'autoFit' not in attrib:
        attrib.update({'autoFit': 'true'})
    
    root = r._root = VisAnt(**attrib)

    t = xml.parse(default_methods_file)
    for elem in t.getroot().findall('method'):
        root.append(method(**elem.attrib))

    root.append(Nodes(**attrib))
    root.append(Edges(**attrib))

    return r

def settag(self, value, name):
    if value != None:
        self.set(name, value)
    else:
        if name in self.attrib:
            del self.attrib[name]

def settext(self, value, name):
    obj = eval(name)
    v = self.find(name)

    if value != None:
        if v:
            self.remove(v)
        self.append(obj(value))
    else:
        if v is not None:
            self.remove(v)

def rgb_list_to_text(lst, method_format=False):
    if method_format:
        sep = ','
    else:
        sep = ' '
    return sep.join([ str(int(x)) for x in lst[:3] ])

def colorwrap(value, method_format=False):
    if isinstance(value, tuple) or isinstance(value, list):
        return rgb_list_to_text(value, method_format)
    if COLORCONVERTER is not None:
        try:
            rgb_tuple = COLORCONVERTER.to_rgb(value)
            rgb_tuple = N.array(rgb_tuple) * 255
            return rgb_list_to_text(rgb_tuple, method_format)
        except:
            pass
    return value

def add_to_set(prop, newelem):
    # list properties return a copy of the internal list, and modification is only supported through direct assignment (i.e., node.groups = grouplist)
    lstcpy = prop
    if lstcpy is not None:
        if newelem in lstcpy:
            return lstcpy

        lstcpy.append(newelem)
    else:
        lstcpy = [newelem]
    return lstcpy

def rem_from_set(prop, name):
    if prop is not None:
        if name in prop:
            lstcpy = set(prop)
            lstcpy.remove(name)
            return list(lstcpy)
    return prop

def get_node_name(node):
    assert isinstance(node, VNodes)
    if node.isduplicate:
        return node.uid
    return node.name


class VisMLTree(object):

    def __init__(self, filename=None):
        self._index = 0

        if filename is not None:
            self.parse(filename)

    def isconnected(self, node1, node2):
        """Takes either node names or node objects as input"""
        return self._edges.isconnected(node1, node2)

    def validate(self):
        for elem in self:
            elem._updatestatic()
            elem.validate()

    def parse(self, filename):
        tree = xml.parse(filename)
        root = tree.getroot()

        assert root.tag == 'VisAnt'

        self._root = VisAnt(**root.attrib)

        for child in root.getchildren():
            self._addbranch(self._root, child)

        for elem in self:
            elem._updatestatic()

    def write(self, filename, prettyprint=True):
        """
        
        Setting the "prettyprint" keyword enables nicely-formatted output via the python xml minidom module

        This is very memory-intensive, so if your xml tree is large, set it to False

        """

        for elem in self:
            elem._updatestatic()
        
        stringout = StringIO.StringIO()
        xml.ElementTree(self._root).write(stringout, xml_declaration=True, encoding="utf-8")
        
        # Hack to allow prettyprint through xml.dom.minidom

        if prettyprint:
            output = minidom.parseString(stringout.getvalue()).toprettyxml()
        else:
            output = stringout.getvalue()

        stringout.close()
        
        f = open(filename, 'w')
        f.write(output)
        f.close()

    @property
    def nodes(self):
        """A list of attached nodes"""
        return self._root.nodes

    @property
    def metanodes(self):
        """A list of attached nodes"""
        return [ x for x in self.nodes if x.ismetanode ]

    @property
    def edges(self):
        """A list of edges in the network"""
        return self._root.edges

    @property
    def _nodes(self):
        """The object itself, rather than a list of nodes"""
        return self._root.find('Nodes')

    @property
    def _edges(self):
        """The object itself, rather than a list of edges"""
        return self._root.find('Edges')

    def __iter__(self):
        return self._root.iter()

    def add_edge(self, node1, node2, method='M0099', **attrib):
        """
        Create a link between node1 and node2

        Note that the link is bidirectional. A VEdge is automatically created.

        Method is required. Default is M0099, which corresponds to 'unknown'.

        """

        try:
            assert isinstance(node1, VNodes) and isinstance(node2, VNodes)
        except:
            raise ValueError, 'Arguments should be of type VNodes'

        node1_name = get_node_name(node1)
        node2_name = get_node_name(node2)

        if not self.isconnected(node1_name, node2_name):
            # Create VEdge
            edge = self._edges.add_edge(node1, node2, **attrib)

        # Add data link
        if node1.isduplicate and node2.isduplicate:
            node1.parent.add_link(node2.parent.name, method, fromDup=node1.uid, toDup=node2.uid, **attrib)
            node2.parent.add_link(node1.parent.name, method, fromDup=node2.uid, toDup=node1.uid, **attrib)
        elif node1.isduplicate or node2.isduplicate:
            dup, nde = node1, node2
            if node2.isduplicate:
                nde, dup = node1, node2

            dup.parent.add_link(nde.name, method, fromDup=dup.uid, **attrib)
            nde.add_link(dup.parent.name, method, toDup=dup.uid, **attrib)
        else:
            node1.add_link(node2.name, method, **attrib)
            node2.add_link(node1.name, method, **attrib)

    def remove_edge(self, node1, node2):
        # TODO This doesn't support unidirectional linkage removal, since add_edge doesn't support adding them
        try:
            assert isinstance(node1, VNodes) and isinstance(node2, VNodes)
        except:
            raise ValueError, 'Arguments should be of type VNodes'

        node1_name = get_node_name(node1)
        node2_name = get_node_name(node2)

        if not self.isconnected(node1_name, node2_name):
            return

        # XXX If there can be multiple edge links via different methods, this function will not clean up properly!

        self._edges.remove_edge(node1, node2)

        if node1.isduplicate:
            node1 = node1.parent
        if node2.isduplicate:
            node2 = node2.parent

        node1.remove_link(node2_name)
        node2.remove_link(node1_name)

    def remove_node_from_metanode(self, node, metanode):
        """
        Removes node from group metanode
        
        Note that in the case where node is a metanode and contains children, these children
        will NOT be added to the parent metanode.

        """

        assert metanode.ismetanode
        node_name = get_node_name(node)

        if node.isduplicate:
            if node.group == metanode.name:
                node.group = None
        else:
            node.remove_group(metanode.name)

        metanode.remove_child(node_name)

    def remove_node(self, node):
        """
        
        Removes the given node.

        Because the node can be
            1) Connected to other nodes
            2) A child of some metanode
            3) A metanode with children
            4) A duplicate
            5) A duplicate which is also the child of a metanode

        Extensive cleanup must be performed

        """
        node_name = get_node_name(node)

        if node.isduplicate:
            node.parent.remove(node)
        else:
            # Delete all duplicates
            if node.duplicates is not None:
                for n in node.duplicates:
                    self.remove_node(n)
        
            # If the node is a metanode, remove all group references from other nodes
            if node.ismetanode:
                for n in self.nodes:
                    n.remove_group(node.name)

            self._nodes.remove(node)

        # Remove all edges connected to the node
        for n in self.nodes:
            if self.isconnected(n, node):
                self.remove_edge(n, node)
    
        # Finally, if any metanodes have this as a child, remove those references
        for n in self.nodes:
            if n.ismetanode:
                n.remove_child(node_name)

    def add_node(self, name, x, y, w='16', h='16', vlabel=None, **attrib):

        index = self._inc_index()

        node = self._nodes.add_node(name, x, y, '0', index, w=w, h=h, vlabel=vlabel, **attrib)
        return node

    def duplicate_node(self, node, x, y, vlabel=None, **attrib):

        index = self._inc_index()
        return node.add_duplicate(x, y, index, vlabel=vlabel, **attrib)

    def add_metanode(self, name, x, y, vlabel=None, children=None, **attrib):

        index = self._inc_index()

        node = self._nodes.add_node(name, x, y, '3', index, vlabel=vlabel, children=children, **attrib)
        return node

    def add_node_to_metanode(self, node, metanode, duplicate=True):
        """
        
        Adds a given node to a given metanode
        
        If 'duplicate' is set to True and the node is already part of a group, 
        the node is duplicated and the duplicate is assigned to the metanode (default VisANT behavior)
        This is ignored if the node is itself a duplicate, and the node will be assigned as normal.

        Do not set duplicate to False unless you know what you are doing!

        As of Feb 27, 2013, there are several VisANT bugs associated with this. Notably, that
        hiding/unhiding a metanode which contains non-duplicated nodes will delete
        said nodes without warning.
        
        """

        assert metanode.ismetanode

        if node.ismetanode:
            metanode.add_child(node.name)
            node.groups = [metanode.name]
            return

        if duplicate and not node.isduplicate and node.groups is not None:
            node = self.duplicate_node(node, node.x, node.y)
        
        node_name = get_node_name(node)
        if node.isduplicate:
            node.group = metanode.name
            node.parent.add_group(metanode.name)
        else:
            node.add_group(metanode.name)

        metanode.add_child(node_name)

        # TODO add_method, add_pathway, add_id, transparent type mechanism (accept/return types, store string internally)
        # Make sure deleting works, plus accounting

    def add_method(self, name, desc, type, visible=None, weight=None, species=None, color=None):
        """

        name: the method id used in VisANT, required. If you need to have your own method, please put the id in the range M7000-M7999 so that VisANT will not enable querying all interaction of this method in the databases.
        desc: the description of the method, required.
        type: the type of the method, either "E" or "C" to indicate experimental or computational, required.
        visible: indicates whether the edge determined by this method only is visible or not, default is true, optional.
        weight: the reliability score of the method, optional, not used so far. 

        """

        self._root.append(method(name=name, desc=desc, type=type, visible=visible, weight=weight, species=species, color=color))

    def _inc_index(self):

        self._index += 1
        return str(self._index)

    def _addbranch(self, root, child):

        if 'index' in child.attrib:
            self._index = max(self._index, int(child.get('index')))

        try:
            obj = eval(child.tag)
        except:
            raise ValueError('Unrecognized tag %s' % child.tag)

        if child.tag in ('method', 'pathway', 'VScreen', 'exLink', 'Nodes', 'VNodes', 'id', 'group', 'Dup', 'link', 'data', 'Edges', 'VEdge'):
            v = obj(**child.attrib)
        elif child.tag in ('vlabel', 'children', 'desc'):
            v = obj(child.text)
        else:
            raise NotImplementedError('Unimplemented tag %s' % child.tag)

        if child.tag == 'Dup':
            v.parent = root

        root.append(v)
        
        for c in child.getchildren():
            self._addbranch(v, c)


class VisMLElement(xml.Element):
    """
    Base class for VisML elements
    Inherits xml.Element
    
    """

    def __init__(self, elementname):

        xml.Element.__init__(self, elementname)

    def __delitem__(self, x):

        warnings.warn('WARNING: Deletion of child elements bugged, try parentelement.remove(child) instead') # TODO

        try:
            self.remove(x)
        except:
            # As long as settag is called, this will result in deleting the attribute
            self.x = None

    def validate(self):
        return True

    def _updatestatic(self):
        pass


class VisAnt(VisMLElement):
    """

    root element

    ver: indicates the version of visML, this attribute is required. However, it is used only for internal control. For most users, simply put an number that is bigger than 1.32 will be fine.
    Note that for compatibility reasons VisAnt will ignore VisML trees that do not begin with ver, e.g., <VisAnt ver=

    species: indicates the species of the network, this attribute is required. If the species of your network is not species-specific, 
             or not in the complete list of species supported in VisANT, simply put "uno" to indicate that the species is "unkown".
             This attribute is useful if you need to query database for further information of your network. For VisANT 2.46 and above, all corresponding database queries will be disabled if the species is "unknown".
    nodecount: indicates the total number of nodes, this attribute is required. Be aware it is the number of total visible and invisible nodes.
    edgeopp: used for meta graph only, optional, default value is false.
    fineArt: used to indicate whether to use better graph drawing, such as thick line, optional, default value is true. Turn it off for large network
    autoFit: used to indicate whether need to fit the network to the default network window, optional, default false.
    double_click: used to change the default behavior when user double-clicking a non-meta-node, if the value of this attribute matches the id attribute of <exLink> double-clicking on the node will launch the link in default browser.
    db_online: used to indicate whether need to query database, optional, default is true. This option will be overwritten if the species is unknown. Available after visML 1.36.
    layout: used to layout the network with specified type of layout, optional. This attribute supports all type of layout that is available in VisANT.
            Available after visML 1.36. Note that this is a read-only attribute. Therefore if a file with this attribute is loaded into VisANT and is later saved by the user,
            this attribute will be gone in the new-saved file. Here is the list of possible values of this attribute: circle, spoke, scramble, relax, embedded, and elegant.
            The later three can have additional parameter for the number of iterations of the layout in the format of embedded:100, as an example. The default iteration is 250.
            Be aware that this number should be bigger if you have a large network. Here is a list of examples with different variants of this optional attribute:
    bgcolor: saved background color
    net_size: Seems to affect the width,height of the network. Changes out the network is displayed, similar to VScreen.

    """

    def __init__(self, ver=None, species=None, net_size=None, bgcolor=None, nodecount=None, edgeopp=None, fineArt=None, autoFit=None, double_click=None, db_online=None, layout=None, **kwds):

        VisMLElement.__init__(self, 'VisAnt')

        # XXX This is a static property now
        #self.nodecount = nodecount

        self.ver       = ver
        self.species   = species
        self.edgeopp   = edgeopp # Bool
        self.fineArt   = fineArt # Bool
        self.autoFit   = autoFit # Bool
        self.layout    = layout  # circle, spoke, scramble, relax, embedded, elegant
        self.double_click = double_click
        self.db_online    = db_online
        self.net_size  = net_size
        self.bgcolor   = bgcolor

    def validate(self):

        try:
            assert all(self.ver, self.species)
        except:
            warn('VisAnt root element is missing required tags (ver, species)')

        if self.double_click:
            try:
                assert self.double_click in [ x.id for x in self.findall('exLink') ]
            except:
                warn('double_click value not found in available exLink elements')

        return True

    @property
    def nodes(self):
        if self.find('Nodes') is not None:
            return self.find('Nodes').findall('VNodes')

    @property
    def edges(self):
        if self.find('Edges') is not None:
            return self.find('Edges').findall('VEdge')

    @property
    def ver(self):
        return self.get('ver')

    @ver.setter
    def ver(self, value):
        try:
            assert float(value) >= 1.32
        except:
            raise ValueError, 'Version should be greater than 1.32; PyVisML does not support earlier versions'

        settag(self, value, 'ver')
    
    @property
    def species(self):
        return self.get('species')

    @species.setter
    def species(self, value):
        # TODO: Assert species is in default list, otherwise 'uno'

        settag(self, value, 'species')

    @property
    def net_size(self):
        return self.get('net_size')

    @net_size.setter
    def net_size(self, value):
        settag(self, value, 'net_size')

    @property
    def bgcolor(self):
        return self.get('bgcolor')

    @bgcolor.setter
    def bgcolor(self, value):
        settag(self, colorwrap(value), 'bgcolor')

    @property
    def nodecount(self):
        if self.nodes is not None:
            return len(self.nodes)
        return 0

    def _updatestatic(self):
        settag(self, str(self.nodecount), 'nodecount')

    @property
    def edgeopp(self):
        return self.get('edgeopp')

    @edgeopp.setter
    def edgeopp(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('edgeopp should be boolean')

        settag(self, value, 'edgeopp')

    @property
    def fineArt(self):
        return self.get('fineArt')

    @fineArt.setter
    def fineArt(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('fineArt should be boolean')

        settag(self, value, 'fineArt')

    @property
    def autoFit(self):
        return self.get('autoFit')

    @autoFit.setter
    def autoFit(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('autoFit should be boolean')

        settag(self, value, 'autoFit')

    @property
    def layout(self):
        return self.get('layout')

    @layout.setter
    def layout(self, value):
        # circle, spoke, scramble, relax, embedded, elegant
        if value is not None:
            v = value.split(':')
            try:
                assert v[0] in ('circle', 'spoke', 'scramble', 'relax', 'embedded', 'elegant')
            except:
                warn('Layout must be of the following types: circle, spoke, scramble, relax, embedded, or elegant.\nThe last three may have iterations specified, e.g., "relax:100"')

        settag(self, value, 'layout')

    @property
    def double_click(self):
        return self.get('double_click')

    @double_click.setter
    def double_click(self, value):
        settag(self, value, 'double_click')

    @property
    def db_online(self):
        return self.get('db_online')

    @db_online.setter
    def db_online(self, value):
        settag(self, value, 'db_online')
        

class method(VisMLElement):
    """

    Method class

    name: the method id used in VisANT, required. If you need to have your own method, please put the id in the range M7000-M7999 so that VisANT will not enable querying all interaction of this method in the databases.
    desc: the description of the method, required.
    type: the type of the method, either "E" or "C" to indicate experimental or computational, required.
    visible: indicates whether the edge determined by this method only is visible or not, default is true, optional.
    weight: the reliability score of the method, optional, not used so far. 

    """

    def __init__(self, name=None, desc=None, type=None, visible=None, weight=None, species=None, color=None):

        VisMLElement.__init__(self, 'method')

        self.name = name
        self.desc = desc
        self.type = type
        self.color   = color
        self.weight  = weight
        self.visible = visible
        self.species = species

    def validate(self):
        try:
            assert all(self.name, self.desc, self.type)
        except:
            warn('Missing required method data (name, desc, type)')

        return True

    @property
    def name(self):
        return self.get('name')
    
    @name.setter
    def name(self, value):
        if value:
            if (value[0] != 'M'):
                warn('Method name should be of the form M####, where # are integers')

            try:
                int(value[1:])
            except:
                warn('Method name should be of the form M####, where # are integers')

        settag(self, value, 'name')

    @property
    def desc(self):
        return self.get('desc')

    @desc.setter
    def desc(self, value):
        settag(self, value, 'desc')

    @property
    def type(self):
        return self.get('type')

    @type.setter
    def type(self, value):
        if value and value not in ('E', 'C'):
            warn('Method type should be E or C')
        settag(self, value, 'type')

    @property
    def visible(self):
        return self.get('visible')

    @visible.setter
    def visible(self, value):
        settag(self, value, 'visible')

    @property
    def weight(self):
        return self.get('weight')

    @weight.setter
    def weight(self, value):
        settag(self, value, 'weight')

    @property
    def species(self):
        return self.get('species')

    @species.setter
    def species(self, value):
        settag(self, value, 'species')

    @property
    def color(self):
        return self.get('color')

    @color.setter
    def color(self, value):
        settag(self, colorwrap(value, method_format=True), 'color')


class pathway(VisMLElement):
    """
    
    Denotes a pathway annotation. New in VisML 1.35

    Has two tags, name (internal ID) and title (a description)
    Both are required

    """

    def __init__(self, name=None, title=None):

        VisMLElement.__init__(self, 'pathway')

        self.name  = name
        self.title = title

    def validate(self):
        try:
            assert all(self.name, self.title)
        except:
            warn('Missing required pathway data (name, title)')

        return True

    @property
    def name(self):
        return self.get('name')
    
    @name.setter
    def name(self, value):
        settag(self, value, 'name')

    @property
    def title(self):
        return self.get('title')
    
    @title.setter
    def title(self, value):
        settag(self, value, 'title')


class VScreen(VisMLElement):
    """
    
    The <VScreen> is an optional element used to retain the zooming level of the network.
    It is not suggest for user to customize it. 
    In case you do need this element, it is highly suggest you use VisANT to zoom at the level you prefer and then save the network.
    You can then copy the elements into your own VisML file.

    x1 y1 x2 y2 w h ps

    """

    def __init__(self, x1=None, y1=None, x2=None, y2=None, w=None, h=None, ps=None):

        VisMLElement.__init__(self, 'VScreen')

        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.w  = w
        self.h  = h
        self.ps = ps

    def validate(self):
        try:
            assert all(self.x1, self.x2, self.y1, self.y2, self.w, self.h, self.ps)
        except:
            warn('Missing required VScreen data (x1, x2, y1, y2, w, h, ps)')

        return True

    @property
    def x1(self):
        return self.get('x1')
    
    @x1.setter
    def x1(self, value):
        settag(self, value, 'x1')

    @property
    def x2(self):
        return self.get('x2')
    
    @x2.setter
    def x2(self, value):
        settag(self, value, 'x2')

    @property
    def y1(self):
        return self.get('y1')
    
    @y1.setter
    def y1(self, value):
        settag(self, value, 'y1')
    
    @property
    def y2(self):
        return self.get('y2')
    
    @y2.setter
    def y2(self, value):
        settag(self, value, 'y2')
    
    @property
    def w(self):
        return self.get('w')
    
    @w.setter
    def w(self, value):
        settag(self, value, 'w')
    
    @property
    def h(self):
        return self.get('h')
    
    @h.setter
    def h(self, value):
        settag(self, value, 'h')
    
    @property
    def ps(self):
        return self.get('ps')
    
    @ps.setter
    def ps(self, value):
        settag(self, value, 'ps')


class exLink(VisMLElement):
    """

    The <exLink> element allows you to add links to external databases for both node and edge, in associated with element <id>.
    When the attribute name of element <id> under <data> matches the attribute name of element <exLink>, a menu will be created in VisANT
    with the name determined by the attribute menu_name in element <exLink> and clicking the menu will launch the default browser
    with the URL determined by the URL attribute of <exLink> and the value attribute of <id> element under <data> element.

    """

    def __init__(self, id=None, menu_name=None, URL=None):

        VisMLElement.__init__(self, 'exLink')

        self.id = None
        self.menu_name = None
        self.URL = None

    def validate(self):
        try:
            assert all(self.id, self.menu_name, self.URL)
        except:
            warn('Missing required exLink data (id, menu_name, URL)')

        # TODO: Validate that there exists a VNodes element with the correct id

        return True

    @property
    def id(self):
        return self.get('id')
    
    @id.setter
    def id(self, value):
        settag(self, value, 'id')

    @property
    def menu_name(self):
        return self.get('menu_name')
    
    @menu_name.setter
    def menu_name(self, value):
        settag(self, value, 'menu_name')

    @property
    def URL(self):
        return self.get('URL')
    
    @URL.setter
    def URL(self, value):
        settag(self, value, 'URL')


class Nodes(VisMLElement):
    """

    Container for VNodes elements

    has a single optional element, size (size of nodes)

    """

    def __init__(self, size=None, **kwds):

        VisMLElement.__init__(self, 'Nodes')

        self.size = None

    @property
    def size(self):
        return self.get('size')
    
    @size.setter
    def size(self, value):
        settag(self, value, 'size')
    
    def add_node(self, name, x, y, type, index, vlabel=None, children=None, w='16', h='16', **attrib):

        attrib.update({'x': x, 'y': y, 'w': w, 'h': h, 'name': name, 'index': index, 'type': type})
        node = VNodes(**attrib)

        settext(node, vlabel, 'vlabel')
        settext(node, children, 'children')

        node.append(data(**attrib))
        self.append(node)

        return node


class VNodes(VisMLElement):
    """
    
    x: x coordinate, required.
    y: y coordinate, required.
    counter: reference counter, its value equals to the number of links connecting to the node. It is also used to determine whether the node is visible or not, set it to at least 1 if the node should be visible, required.
    w: width of the node, required.
    h: height of the node, required.
    labelOn: determine whether the node label is visible, optional, default false.
    linkDisplayed: used to indicate whether all the nodes link is displayed, special designed for step-wise expansion of the interaction network, optional, default false. When it is true, the node show "-" sign, otherwise, "+" sign.
    linkLoaded: indicated whether the node has been queried against database , optional, default false. Set it to be true if you do not want to query database when double-click the node.
    extend: designed for step-wise expansion of the interaction network, optional, default true. Set it to true if you do not want the node position to be changed when double-click the node.
    size: node size, optional. Default -1:auto, Range: 4-30.
    ncc: RGB value of node color, optional.
    labelSize: size of the node label, optional. Default -1:auto, Range: 6-25.
    labelPos: position of the label, optional. Default -1:auto, 0:center, 1:left, 2:right, 3:above, 4:below
    esymbol: determine whether to show the extension symbol (+/-) or not, optional, default is true.
    ksymbol: 
    group: if the node is in a group, this attribute represents the corresponding group id, optional. 
    exShape: metanode polygon shape, optional. Current options are None: auto, 2001: convex_polygon, and 2002: round_rectangle
    properties not present in the xml:
    visible: changes the counter to negative, which prevents the node from being shown
    
    """

    def __init__(self, x=None, y=None, counter=None, w=None, h=None, exShape=None, labelOn=None, linkDisplayed=None, linkLoaded=None, extend=None, size=None, ncc=None, labelStyle=None, labelSize=None, labelPos=None, esymbol=None, ksymbol=None, group=None, axisnode=None, visible=True, eletag='VNodes', dt=None, childVisible=None, **kwds):

        # We make an exception here so Dup can completely inherit VNodes methods and properties
        VisMLElement.__init__(self, eletag)

        # XXX Static property
        # self.counter = counter

        self.x = x
        self.y = y
        self.w = w
        self.h = h
        self.labelOn       = labelOn       # Bool
        self.linkDisplayed = linkDisplayed # Bool
        self.linkLoaded    = linkLoaded    # Bool
        self.extend        = extend        # Bool
        self.esymbol       = esymbol       # Bool
        self.ksymbol       = ksymbol       # Bool
        self.size = size
        self.ncc  = colorwrap(ncc)
        self.labelSize = labelSize
        self.labelStyle = labelStyle
        self.labelPos  = labelPos
        self.group     = group
        self.axisnode  = axisnode
        self.visible   = visible
        self.exShape   = exShape
        self.duplicate_index = -1
        self.dt   = dt
        self.childVisible = childVisible

    def validate(self):
        try:
            assert all(self.x, self.y, self.w, self.h)
        except:
            warn('Missing required VNodes data (x, y, w, h)')

        # TODO: Validate that there are metanodes for all claimed group ids

        return True

    def add_duplicate(self, x, y, index, vlabel=None, w='16', h='16', **attrib):
        if self.isduplicate:
            raise ValueError, 'Unable to duplicate a duplicate node'
        
        idx = self._inc_dup_index()
        uid = '_'.join([self.name, idx])

        attrib.update({'x': x, 'y': y, 'w': w, 'h': h, 'uid': uid, 'index': index, 'parent': self})

        dup = Dup(**attrib)
        settext(dup, vlabel, 'vlabel')

        self.append(dup)
        return dup

    def add_group(self, name):
        """Adds the specified group name to groups. Does not change the group node children! Use VisMLTree.add_node_to_metanode instead"""
        if self.isduplicate:
            raise ValueError, 'Unable to add group to duplicate node. Use node.group = value instead'
        self.data.add_group(name)

    def remove_group(self, name):
        """Adds the specified group name to groups. Does not change the group node children! Use VisMLTree.remove_node_from_metanode instead"""
        if self.isduplicate:
            raise ValueError, 'Unable to remove group from duplicate node. Use node.group = value instead'
        self.data.remove_group(name)

    def add_pathway(self, name):
        if self.isduplicate:
            raise ValueError, 'Unable to set pathways for a duplicate node'
        self.data.add_pathway(name)

    def remove_pathway(self, name):
        if self.isduplicate:
            raise ValueError, 'Unable to set pathways for a duplicate node'
        self.data.remove_pathway(name)

    def add_child(self, name):
        if self.isduplicate:
            raise ValueError, 'Unable to set children for a duplicate node'
        self.children = add_to_set(self.children, name)

    def remove_child(self, name):
        if self.isduplicate:
            raise ValueError, 'Unable to set children for a duplicate node'
        self.children = rem_from_set(self.children, name)

    def add_link(self, to, method, **attrib):
        self.data.add_link(to, method, **attrib)

    def remove_link(self, name):
        self.data.remove_link(name)

    def _inc_dup_index(self):
        self.duplicate_index += 1
        return str(self.duplicate_index)

    @property
    def isduplicate(self):
        if 'uid' in self.attrib:
            return True
        return False

    @property
    def type(self):
        return self.data.type

    @property
    def links(self):
        return self.data.links

    @property
    def children(self):
        if self.find('children') is not None:
            if self.find('children').text is not None:
                return self.find('children').text.split(',')[:-1] # format is CHILD1,CHILD2,

    @children.setter
    def children(self, value):
        assert isinstance(value, list) or value is None

        if not value:
            if self.children is not None:
                self.remove(self.find('children'))
            return
    
        text = ','.join(value + [''])      # format is CHILD1,CHILD2,

        if self.children is not None:
            self.find('children').text = text
        else:
            self.append(children(text))

    @property
    def data(self):
        if self.isduplicate:
            raise ValueError, 'Duplicate node has no data container'
        return self.find('data')

    @property
    def ismetanode(self):
        return self.data.ismetanode

    @property
    def duplicates(self):
        return self.findall('Dup')

    @property
    def pathways(self):
        return self.data.pathways

    @pathways.setter
    def pathways(self, value):
        self.data.pathways = value

    @property
    def desc(self):
        return self.data.desc

    @desc.setter
    def desc(self, value):
        self.data.desc = value

    @property
    def groups(self):
        if self.isduplicate:
            return [self.group]
        return self.data.groups

    @groups.setter
    def groups(self, value):
        if self.isduplicate:
            if isinstance(value, list):
                raise ValueError, 'Trying to set a group list for a duplicate node!'
            self.group = value
            return
        self.data.groups = value

    @property
    def name(self):
        return self.data.name

    @name.setter
    def name(self, value):
        # TODO: This should regenerate all the links to this node
        self.data.name = value

    @property
    def vlabel(self):
        if self.find('vlabel') is not None:
            return self.find('vlabel').text

    @vlabel.setter
    def vlabel(self, value):
        settext(self, value, 'vlabel')

    @property
    def duplicates(self):
        return self.findall('Dup')

    @property
    def counter(self):
        # TODO: Real data
        if self.visible:
            return 1
        return -1

    def _updatestatic(self):
        settag(self, str(self.counter), 'counter')
        
        # XXX findall has reported [] instead of None in certain test cases
        if self.duplicates:
            self.duplicate_index = max([ int(uid.split('_')[-1]) for uid in [ dup.uid for dup in self.duplicates ] ])

    @property
    def x(self):
        return self.get('x')
    
    @x.setter
    def x(self, value):
        settag(self, value, 'x')

    @property
    def y(self):
        return self.get('y')
    
    @y.setter
    def y(self, value):
        settag(self, value, 'y')

    @property
    def w(self):
        return self.get('w')
    
    @w.setter
    def w(self, value):
        settag(self, value, 'w')

    @property
    def h(self):
        return self.get('h')
    
    @h.setter
    def h(self, value):
        settag(self, value, 'h')
    
    @property
    def childVisible(self):
        return self.get('childVisible')
    
    @childVisible.setter
    def childVisible(self, value):
        settag(self, value, 'childVisible')

    @property
    def axisnode(self):
        return self.get('axisnode')
    
    @axisnode.setter
    def axisnode(self, value):
        settag(self, value, 'axisnode')

    @property
    def labelOn(self):
        return self.get('labelOn')
    
    @labelOn.setter
    def labelOn(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('labelOn should be boolean')

        settag(self, value, 'labelOn')

    @property
    def exShape(self):
        return self.get('exShape')
    
    @exShape.setter
    def exShape(self, value):
        try:
            if value is not None:
                assert value in ['2001', '2002']
        except:
            warn('Current exShape support includes None, "2001", and "2002", which is auto, convex_polygon, and round_rectangle, respectively')
            value = None

        settag(self, value, 'exShape')

    @property
    def linkDisplayed(self):
        return self.get('linkDisplayed')
    
    @linkDisplayed.setter
    def linkDisplayed(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('linkDisplayed should be boolean')

        settag(self, value, 'linkDisplayed')

    @property
    def linkLoaded(self):
        return self.get('linkLoaded')
    
    @linkLoaded.setter
    def linkLoaded(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('linkLoaded should be boolean')

        settag(self, value, 'linkLoaded')

    @property
    def size(self):
        return self.get('size')
    
    @size.setter
    def size(self, value):
        settag(self, value, 'size')

    @property
    def extend(self):
        return self.get('extend')
    
    @extend.setter
    def extend(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('extend should be boolean')

        settag(self, value, 'extend')

    @property
    def esymbol(self):
        return self.get('esymbol')
    
    @esymbol.setter
    def esymbol(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('esymbol should be boolean')

        settag(self, value, 'esymbol')
    
    @property
    def ksymbol(self):
        return self.get('ksymbol')
    
    @ksymbol.setter
    def ksymbol(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('ksymbol should be boolean')

        settag(self, value, 'ksymbol')

    @property
    def ncc(self):
        return self.get('ncc')
    
    @ncc.setter
    def ncc(self, value):
        settag(self, colorwrap(value), 'ncc')

    @property
    def dt(self):
        return self.get('dt')

    @dt.setter
    def dt(self, value):
        # Relevant VisANT codeblock:
        # //node shape
        # public int CIRCLE=0,RECTANGLE=4, VBALL=62, CIRCLE_COMPOUND=63,TRIANGLE=64, DIAMOND=65, CIRCLE_DRUG=69;
        # public int HEXAGON=66, OCTAGON=67, SQUARE=68, ROUND_RECTANGLE=3, RECTANGLE_3D=1;
        # //10-19 resevered for the shape fitted for the label
        # public int RECTANGLE_FIT=10,ROUND_RECTANGLE_FIT=11,RECTANGLE_3D_FIT=19;
        # public int EXP_CONVEXITY=2001, EXP_RECTANGLE=2002, EXP_AUTO=-1;
        # public int EXPRESSION_PLOT=1000;

        shapes = {'circle': 0, 'rectangle': 4, 'vball': 62, 'circle_compound': 63, 'triangle': 64,
                'diamond': 65, 'circle_drug': 69, 'hexagon': 66, 'octagon': 67, 'square': 68,
                'round_rectangle': 3, 'rectangle_3D': 1, 'rectangle_fit': 10, 'round_rectangle_fit': 11,
                'rectangle_3D_fit': 19, 'exp_convexity': 2001, 'exp_rectangle': 2002, 'exp_auto': -1,
                'expression_plot': 1000}

        try:
            v = str(int(value))
        except:
            # Shape name given

            if value is not None:
                try:
                    assert value in shapes
                    v = str(shapes[value])
                except:
                    warn('Unknown shape. Allowed values: ' + ', '.join(shapes.keys()) + '\n')
                    v = None
            else:
                v = None

        settag(self, v, 'dt')

    @property
    def labelPos(self):
        return self.get('labelPos')
    
    @labelPos.setter
    def labelPos(self, value):
        settag(self, value, 'labelPos')

    @property
    def labelSize(self):
        return self.get('labelSize')
    
    @labelSize.setter
    def labelSize(self, value):
        settag(self, value, 'labelSize')

    @property
    def labelStyle(self):
        return self.get('labelStyle')
    
    @labelStyle.setter
    def labelStyle(self, value):
        settag(self, value, 'labelStyle')


    @property
    def group(self):
        return self.get('group')
    
    @group.setter
    def group(self, value):

        # TODO: Validate the existence of a metanode to assign this to?
        settag(self, value, 'group')


class vlabel(VisMLElement):
    """
    Optional node label, no children or tags (only text)
    
    """
    def __init__(self, text):
        VisMLElement.__init__(self, 'vlabel')

        self.text = text

class children(VisMLElement):
    """
    Metanode child container, no children or tags (only text)
    
    """
    def __init__(self, text):
        VisMLElement.__init__(self, 'children')

        self.text = text

class desc(VisMLElement):
    """
    data or link child container, no children or tags (only text)
    
    """
    def __init__(self, text):
        VisMLElement.__init__(self, 'desc')

        self.text = text

class id(VisMLElement):
    """    
    child of data or link

    the <id> element can also be used to create HTTP link by matching to the id attribute of the <exLink> element.

    """
    def __init__(self, name=None, value=None):
        
        VisMLElement.__init__(self, 'id')

        self.name  = name
        self.value = value

    def validate(self):

        try:
            assert all(self.name, self.value)
        except:
            warn('Element id missing required tags (self.name, self.value)')

        return True

    @property
    def name(self):
        return self.get('name')
    
    @name.setter
    def name(self, value):
        settag(self, value, 'name')
        
    @property
    def value(self):
        return self.get('value')
    
    @value.setter
    def value(self, value):
        settag(self, value, 'value')


class group(VisMLElement):
    """
    child of data

    These elements will appear if the node represented by the data element is in a group. This element has two attributes: name determines the type of group ("common" for most cases),
    while the value represents the group ids delimited by ",", in case there are more than one group (in such case, the node itself must have duplications).

    """
    def __init__(self, name=None, value=None):
        
        VisMLElement.__init__(self, 'group')

        self.name  = name
        self.value = value

    def validate(self):

        try:
            assert all(self.name, self.value)
        except:
            warn('Element group missing required tags (self.name, self.value)')

        return True

    @property
    def name(self):
        return self.get('name')
    
    @name.setter
    def name(self, value):
        settag(self, value, 'name')
        
    @property
    def value(self):
        return self.get('value')
    
    @value.setter
    def value(self, value):
        settag(self, value, 'value')


class Dup(VNodes):
    """
    In case of duplication, this element will be appear under <VNodes> element. this element describes visual properties of the duplicated nodes and have similar attributes as <VNodes>, except following:
    uid: the identification of the duplicated node, required.
    group: if the duplicated node is in a group, this attribute denotes the corresponding group id. optional. 

    """
    def __init__(self, uid=None, index=None, x=None, y=None, counter=None, w=None, h=None, labelOn=None, linkDisplayed=None, linkLoaded=None,
                 extend=None, size=None, ncc=None, labelSize=None, labelPos=None, esymbol=None, ksymbol=None, group=None, visible=True, parent=None):
        
        VNodes.__init__(self, x=x, y=y, counter=counter, w=w, h=h, labelOn=labelOn, linkDisplayed=linkDisplayed, linkLoaded=linkLoaded, extend=extend,
                        size=size, ncc=ncc, labelSize=labelSize, labelPos=labelPos, esymbol=esymbol, ksymbol=ksymbol, group=group, visible=visible, eletag='Dup')

        self.parent = parent
        self.uid = uid
        self.index = index

    def validate(self):

        try:
            assert all(self.uid)
        except:
            warn('Dup nodes require uid to be specified')

        return True

    @property
    def uid(self):
        return self.get('uid')
    
    @uid.setter
    def uid(self, value):
        settag(self, value, 'uid')

    @property
    def index(self):
        return self.get('index')
    
    @index.setter
    def index(self, value):
        settag(self, value, 'index')


class link(VisMLElement):
    """
    
    child of data

    This element has following attributes:
    to: the node name that this link connects to, required.
    toDup: If this node points to a specific duplicate in the target, its uid should be here
    fromDup: If a duplicate in this node points to a target, its uid should be here
    method: the method id associated with the link. If there is a literature reference corresponding pubmed id can be appended as the format shown below. required.
    toType: integer number to determine the type of the to-end of the edge, optional, default is 0. Its value can be ranged (please reference KEGG database for the meaning of different edge type:
    fromType: integer number to determine the type of the from-end of the edge, optional. The value has the exact same range as toType.
    weight: value for edge weight. Can specify multiple method weights with weight=MXXXX:value, currently unsupported by PyVisML

    """
    
    # XXX Support method weights

    def __init__(self, to=None, method=None, toType=None, fromType=None, weight=None, desc=None, toDup=None, fromDup=None, **kwds):

        VisMLElement.__init__(self, 'link')

        self.to = to
        self.method = method
        self.toType = toType
        self.fromType = fromType
        self.desc = desc
        self.fromDup = fromDup
        self.toDup = toDup
        self.weight = weight

    def validate(self):

        try:
            assert all(self.to, self.method)
        except:
            warn('link element missing required data (to, method)')

        return True

    @property
    def target(self):
        if self.toDup is not None:
            return self.toDup
        return self.to

    @property
    def desc(self):
        if self.find('desc') is not None:
            return self.find('desc').text

    @desc.setter
    def desc(self, value):
        settext(self, value, 'desc')

    @property
    def to(self):
        return self.get('to')
    
    @to.setter
    def to(self, value):
        settag(self, value, 'to')

    @property
    def method(self):
        return self.get('method')
    
    @method.setter
    def method(self, value):
        if value:
            if (value[0] != 'M'):
                warn('Method name should be of the form M####, where # are integers')

            try:
                int(value[1:])
            except:
                warn('Method name should be of the form M####, where # are integers')

        settag(self, value, 'method')

    @property
    def toType(self):
        return self.get('toType')
    
    @toType.setter
    def toType(self, value):
        settag(self, value, 'toType')
        
    @property
    def fromDup(self):
        return self.get('fromDup')
    
    @fromDup.setter
    def fromDup(self, value):
        settag(self, value, 'fromDup')
    
    @property
    def toDup(self):
        return self.get('toDup')
    
    @toDup.setter
    def toDup(self, value):
        settag(self, value, 'toDup')

    @property
    def fromType(self):
        return self.get('fromType')
    
    @fromType.setter
    def fromType(self, value):
        settag(self, value, 'fromType')

    @property
    def weight(self):
        return self.get('weight')
    
    @weight.setter
    def weight(self, value):
        settag(self, value, 'weight')


class data(VisMLElement):
    """
    VNode data child

    This element has following attributes:
    name: default node name, required.
    index: integer used to identify node, should be different for each <data> element, including <Dup>, required.
    type: type of the node, optional, default 0. built in value=0:normal protein/gene, 1:chemical compound, 2:KEGG Map,3: general group node, 4:Protein Complex, 5:Node of Go Term, 6:Pathway Node, 7:KEGG group node, >=100:user added
    alias: a list of alias delimited ",".

    The <data> element may have following child elements:
    <desc>: description of the node, will be shown in the tooltip when mouse-over the node, optional.
    <id>: external database ids of the node, optional. If its name attribute matches the id attribute of <exLink> element, 
    a menu will be created to allow the http link to external database as mentioned earlier. There can be more than one id elements.
    <link>: this element represents a link starts from the node represented by element, and stops at the node indicated by the to attribute, optional (a node can have no links/edges). 
    It needs to be aware that one edge may have more than one links.
    <group>: this elements will appear if the node represented by the data element is in groups. 
    This element has two attribute: name determines the type of group ("common" for most cases), 
    while the value represents the group ids delimited by ",", in case there are more than one group (in such case, the node itself must have duplications).

    """

    def __init__(self, name=None, index=None, type=None, alias=None, pws=None, **kwds):

        VisMLElement.__init__(self, 'data')

        self.ismetanode = None # Set this first because self.type assignment will overwrite it
        self.name  = name
        self.index = index
        self.type  = type
        self.alias = alias
        self.pws   = pws

    def validate(self):
        try:
            assert all(self.name, self.index)
        except:
            warn('Missing required data tags (name, index)')

        if self.type is not None:
            if self.metanode:
                pass
                # TODO: Validate that it has children

        # TODO: Assert index is unique

        return True
    
    def add_link(self, to, method, **attrib):
        # TODO: dict
        target = to
        if 'toDup' in attrib:
            target = attrib['toDup']

        if self.links is not None:
            for edge in self.links:
                if edge.target == target:
                    return

        self.append(link(to=to, method=method, **attrib))

    def remove_link(self, target):
        # TODO: dict
        if self.links is not None:
            for edge in self.links:
                if edge.target == target:
                    self.remove(edge)
                    return

    def remove_group(self, name):
        self.groups = rem_from_set(self.groups, name)
    
    def remove_pathway(self, name):
        self.pathways = rem_from_set(self.pathways, name)

    def add_group(self, name):
        self.groups = add_to_set(self.groups, name)

    def add_pathway(self, name):
        self.pathways = add_to_set(self.pathways, name)

    @property
    def _group_dict(self):
        """Returns a dict of group children, where keys are 'name' tags and values are the group children elements themselves"""
        if self.find('group') is not None:
            return dict([ (grp.name, grp) for grp in self.findall('group') ])

    @property
    def desc(self):
        if self.find('desc') is not None:
            return self.find('desc').text

    @desc.setter
    def desc(self, value):
        settext(self, value, 'desc')

    @property
    def aliases(self):
        if self.alias:
            return self.alias.split(',')

    def _set_groups(self, name, value):
        if value:
            if self._group_dict is not None and name in self._group_dict:
                self._group_dict[name].value = value
            else:
                self.append(group(name=name, value=value))
        else:
            if self._group_dict is not None and name in self._group_dict:
                self.remove(self._group_dict[name])

    @property
    def pathways(self):
        """A list of pathways the node belongs to. Note that the data tag 'pws' will override any 'group' child element of data which has 'pathways' as its name"""
        if self.pws:
            return self.pws.split(' ')

    @pathways.setter
    def pathways(self, value):
        
        assert isinstance(value, list) or value is None
        plist = value
        if value is not None:
            plist = ' '.join(value)

        settag(self, plist, 'pws')
        self._set_groups('pathway', plist)

    @property
    def groups(self):
        """A list of groups the node belongs to, if any."""
        if self._group_dict and 'common' in self._group_dict:
            return self._group_dict['common'].value.split(',')

    @groups.setter
    def groups(self, value):
        
        assert isinstance(value, list) or value is None
        glist = value
        if value is not None:
            glist = ','.join(value)
        
        self._set_groups('common', glist)

    @property
    def links(self):
        return self.findall('link')

    @property
    def ids(self):
        return self.findall('id')

    @property
    def type(self):
        return self.get('type')
    
    @type.setter
    def type(self, value):
        try:
            v = int(value)
        except:
            warn('Type should be an integer')

        if v in range(3, 8):
            self.ismetanode = True
        settag(self, value, 'type')

    @property
    def pws(self):
        return self.get('pws')

    @pws.setter
    def pws(self, value):
        settag(self, value, 'pws')

    @property
    def name(self):
        return self.get('name')
    
    @name.setter
    def name(self, value):
        #if isinstance(value, str):
        #    value = value.upper()
        settag(self, value, 'name')

    @property
    def index(self):
        return self.get('index')
    
    @index.setter
    def index(self, value):
        settag(self, value, 'index')

    @property
    def alias(self):
        return self.get('alias')
    
    @alias.setter
    def alias(self, value):
        settag(self, value, 'alias')


class Edges(VisMLElement):
    """
    
    edges are listed under <Edges> element named as VEdge.     

    Element <Edges> has following attributes:
    thick_line: the option to use the line thickness to represent the edge weight, default is true, optional.
    color_line: the option to use the line color to represent the edge weight, default is false, optional.
    exclude_meta: the option to exclude the metaedges when use weight cutoff, optional.

    Following lists an example of this element with all attributes on.
    <Edges thick_line="false" color_line="true" exclude_meta="true">

    """

    def __init__(self, thick_line=None, color_line=None, exclude_meta=None, opacity=None, **kwds):

        VisMLElement.__init__(self, 'Edges')
            
        self.thick_line = thick_line
        self.color_line = color_line
        self.exclude_meta = exclude_meta
        self.opacity = opacity
        self._edge_dict = {}

    def validate(self):

        # TODO: assert there are edges?
        return True

    def isconnected(self, node1, node2):
        """Returns true if a link between node names name1 and name2 is found in the edge list"""

        node1_name, node2_name = node1, node2
        
        if isinstance(node1, VNodes):
            node1_name = get_node_name(node1)
        if isinstance(node2, VNodes):
            node2_name = get_node_name(node2)

        if node1_name in self._edge_dict and node2_name in self._edge_dict[node1_name]:
            return True
        return False

    def add_edge(self, node1, node2, **attrib):
        """
        Helper function to create an edge between node1 and node2
        Updates the internal edge dictionary
        
        """

        if isinstance(node1, VNodes):
            node1_name = get_node_name(node1)
        if isinstance(node2, VNodes):
            node2_name = get_node_name(node2)

        if not self.isconnected(node1_name, node2_name):
            attrib.update({'linkFrom': node1_name, 'to': node2_name})
            edge = VEdge(**attrib)
            self.append(edge)
            
            self._edge_dict.setdefault(node1_name, set()).add(node2_name)
            self._edge_dict.setdefault(node2_name, set()).add(node1_name)

    def remove_edge(self, node1, node2):
        """
        Helper function to remove an edge between node1 and node2
        Updates the internal edge dictionary

        """

        if isinstance(node1, VNodes):
            node1_name = get_node_name(node1)
        if isinstance(node2, VNodes):
            node2_name = get_node_name(node2)

        if self.isconnected(node1_name, node2_name):
            for edge in self.edges:
                if (edge.linkFrom == node1_name and edge.to == node2_name) or (edge.linkFrom == node2_name and edge.to == node1_name):
                    self.remove(edge)
                    self._edge_dict[node1_name].remove(node2_name)
                    self._edge_dict[node2_name].remove(node1_name)
                    return

    @property
    def edges(self):
        return self.findall('VEdge')

    @property
    def thick_line(self):
        return self.get('thick_line')
    
    @thick_line.setter
    def thick_line(self, value):
        settag(self, value, 'thick_line')

    @property
    def opacity(self):
        return self.get('opacity')

    @opacity.setter
    def opacity(self, value):
        settag(self, value, 'opacity')

    @property
    def color_line(self):
        return self.get('color_line')
    
    @color_line.setter
    def color_line(self, value):
        settag(self, value, 'color_line')

    @property
    def exclude_meta(self):
        return self.get('exclude_meta')
    
    @exclude_meta.setter
    def exclude_meta(self, value):
        settag(self, value, 'exclude_meta')

    def _updatestatic(self):
        ed  = {}

        for edge in self.edges:
            ed.setdefault(edge.linkFrom, set()).add(edge.to)
            ed.setdefault(edge.to, set()).add(edge.linkFrom)

        self._edge_dict = ed


class VEdge(VisMLElement):
    """
    
    edges are listed under <Edges> element named as VEdge.     

    Element <VEdge> has following attributes:

    from: the id of from node, required. Renamed to linkFrom to stop collisions with python keywords
    to: the id of to node, required.
    elabel: the label of the edge, optional.
    la: boolean flag to determine whether the edge label shown be visible, optional, default false. 

    """

    def __init__(self, linkFrom=None, to=None, elabel=None, la=None, **kwds):

        VisMLElement.__init__(self, 'VEdge')
        
        if 'from' in kwds:
            self.linkFrom = kwds['from']
        else:
            self.linkFrom = linkFrom

        self.to       = to
        self.elabel   = elabel
        self.la       = la

    def validate(self):
        try:
            assert all(self.linkFrom, self.to)
        except:
            warn('Element VEdge missing required tags (linkFrom [from], to)')

        # TODO: Validate link from and to nodes as existent

        return True

    @property
    def linkFrom(self):
        return self.get('from')
    
    @linkFrom.setter
    def linkFrom(self, value):
        settag(self, value, 'from')

    @property
    def to(self):
        return self.get('to')
    
    @to.setter
    def to(self, value):
        settag(self, value, 'to')

    @property
    def elabel(self):
        return self.get('elabel')
    
    @elabel.setter
    def elabel(self, value):
        settag(self, value, 'elabel')

    @property
    def la(self):
        return self.get('la')
    
    @la.setter
    def la(self, value):
        try:
            if value is not None:
                assert value in ['true', 'false']
        except:
            warn('la should be boolean')
        settag(self, value, 'la')

