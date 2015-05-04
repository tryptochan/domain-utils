"""Module to work with network in SIF format.

http://wiki.cytoscape.org/Cytoscape_User_Manual/Network_Formats#SIF_Format

"""

class SIFNode(object):
    def __init__(self, id):
        self.id = id
        self.edges = []
        self.reverse= []


class SIFEdge(object):
    def __init__(self, type, target):
        self.type = type
        self.target = target


class SIFNetwork(object):

    def __init__(self):
        self.nodes = {}

    def parse(self, fp):
        if isinstance(fp, basestring):
            fp = open(fp, 'r')
        for l in fp:
            tmp = l.split()
            source = tmp[0]
            edge_type = tmp[1]
            if source not in self.nodes:
                source_node = SIFNode(source)
                self.nodes[source] = source_node
            else:
                source_node = self.nodes[source]
            for target in tmp[2:]:
                if target not in self.nodes:
                    target_node = SIFNode(target)
                    self.nodes[target] = target_node
                else:
                    target_node = self.nodes[target]
                edge = SIFEdge(edge_type, target_node)
                source_node.edges.append(edge)
                r_edge = SIFEdge(edge_type, source_node)
                target_node.reverse.append(r_edge)
        fp.close()

    def subnetwork(self, center_id, degree=1):
        """Find the subnetwork with up to a given degree neighbors,
        regardless of direction.
        """
        sub = self._subnetwork(center_id, degree)
        for n in sub.nodes.itervalues():
            for edge in n.edges:
                complement_flag = False
                for r_edge in edge.target.reverse:
                    if n.id == r_edge.target.id and edge.type == r_edge.type:
                        complement_flag = True
                        break
                if not complement_flag:
                    edge.target.reverse.append(SIFEdge(edge.type, n))
            for edge in n.reverse:
                complement_flag = False
                for r_edge in edge.target.edges:
                    if n.id == r_edge.target.id and edge.type == r_edge.type:
                        complement_flag = True
                        break
                if not complement_flag:
                    edge.target.edges.append(SIFEdge(edge.type, n))

        return sub

    def _subnetwork(self, id, degree=1, sub=None, visited=None):
        if degree == 0:
            return None
        if visited is None:
            visited = set()
        if id in visited:
            return None
        if sub is None:
            sub = SIFNetwork()
        if id not in sub.nodes:
            center = SIFNode(id)
            sub.nodes[id] = center
        else:
            center = sub.nodes[id]

        for edge in self.nodes[id].edges:
            target_node = edge.target
            if target_node.id in visited:
                continue
            if target_node.id not in sub.nodes:
                new_target = SIFNode(target_node.id)
                sub.nodes[target_node.id] = new_target
            else:
                new_target = sub.nodes[target_node.id]
            center.edges.append(SIFEdge(edge.type, new_target))
            self._subnetwork(target_node.id, degree-1, sub, visited)

        for edge in self.nodes[id].reverse:
            target_node = edge.target
            if target_node.id in visited:
                continue
            if target_node.id not in sub.nodes:
                new_target = SIFNode(target_node.id)
                sub.nodes[target_node.id] = new_target
            else:
                new_target = sub.nodes[target_node.id]
            center.reverse.append(SIFEdge(edge.type, new_target))
            self._subnetwork(target_node.id, degree-1, sub, visited)

        visited.add(id)
        return sub

    def write(self, filename):
        with open(filename, 'w') as fp:
            for node in self.nodes.itervalues():
                for edge in node.edges:
                    fp.write('%s\t%s\t%s\n' %
                            (node.id, edge.type, edge.target.id))

    def to_json(self):
        """For use with Cytoscape.js"""
        import json
        js = []
        for node in self.nodes.itervalues():
            js.append({'group': 'nodes',
                      'data': {'id': node.id}
                      })
            for edge in node.edges:
                js.append({'group': 'edges',
                          'data': {
                              'source': node.id,
                              'target': edge.target.id,
                              'type': edge.type
                              }
                          })
        return json.dumps(js)


# vim: ts=4 expandtab sw=4 sts=4 tw=78
