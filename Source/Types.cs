namespace LeidenCD;

using System.Collections.Generic;
using System.Diagnostics.CodeAnalysis;
using System.Linq;

public readonly record struct NodeId(int Id);
public readonly record struct EdgeId(int Id);

public sealed partial class Graph<N, E>
{
    readonly List<N> _Nodes = [];
    readonly List<EdgeInfo> _Edges = [];
    internal readonly List<Dictionary<int, int>> _Connections = [];
    float _AllWeight = 0.0f;

    public NodeId AddNode(N node)
    {
        var id = new NodeId(_Nodes.Count);
        _Nodes.Add(node);
        _Connections.Add([]);
        return id;
    }

    public EdgeId? AddEdge(NodeId node1, NodeId node2, E edge, float weight)
    {
        if (node1.Id == node2.Id)
        {
            return null;
        }

        // normalize the connection lookup keys
        if (node1.Id > node2.Id)
        {
            (node1, node2) = (node2, node1);
        }

        // check if the edge already exists
        var conn = _Connections[node1.Id];
        if (conn.TryGetValue(node2.Id, out var edgeId))
        {
            var old = _Edges[edgeId];
            _Edges[edgeId] = old with { Weight = weight + old.Weight, Data = edge };
            _AllWeight += weight;
            return new EdgeId(edgeId);
        }

        edgeId = conn[node2.Id] = _Edges.Count;
        var edgeInfo = new EdgeInfo { Weight = weight, Data = edge, EdgeId = edgeId };
        _Edges.Add(edgeInfo);
        _AllWeight += weight;

        return new EdgeId(edgeId);
    }

    public NodeEdges GetNodeEdges(NodeId id) => new(this, id);

    public record struct EdgeInfo
    {
        public required float Weight;
        public required E Data;
        public required int EdgeId;
    }

    public int CountNodes() => _Nodes.Count;

    public N GetNodeData(NodeId id) => _Nodes[id.Id];
    public EdgeInfo GetEdgeInfo(EdgeId id) => _Edges[id.Id];
    public EdgeInfo? GetEdgeInfo(NodeId node1, NodeId node2)
    {
        // normalize the connection lookup keys
        if (node1.Id > node2.Id)
        {
            (node1, node2) = (node2, node1);
        }
        return _Connections[node1.Id].TryGetValue(node2.Id, out var edgeId) ? _Edges[edgeId] : null;
    }

    public readonly struct NodeEdges(Graph<N, E> graph, NodeId id)
    {
        readonly Graph<N, E> _Graph = graph;
        readonly NodeId _Id = id;

        public IEnumerable<EdgeId> IterEdges() => _Graph._Connections[_Id.Id].Select(kvp => new EdgeId(kvp.Value));
        public int CountEdges() => _Graph._Connections[_Id.Id].Count;
        public bool TryGetEdge(NodeId node, [MaybeNullWhen(false)] out EdgeId? edge)
        {
            if (_Graph._Connections[_Id.Id].TryGetValue(node.Id, out var edgeId))
            {
                edge = new EdgeId(edgeId);
                return true;
            }
            edge = default;
            return false;
        }
    }
}
