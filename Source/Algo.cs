using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Threading.Tasks;

namespace LeidenCD;

/*
    type Community =
        | L1Community of HashSet<int>
        | LNCommunity of Community[]
*/
public abstract record Community
{
    public int CountNodes()
    {
        Stack<Community> stack = new();
        stack.Push(this);
        int count = 0;
        while (stack.TryPop(out var community))
        {
            switch (community)
            {
                case L1Community l1:
                    count += l1.Nodes.Count;
                    break;
                case LNCommunity ln:
                    foreach (var subCommunity in ln.Communities)
                    {
                        stack.Push(subCommunity);
                    }
                    break;
                default:
                    throw new InvalidOperationException($"Unknown community type: {community.GetType().Name}");
            }
        }

        return count;
    }


    public IEnumerable<int> FetchAllNodes()
    {
        if (this is L1Community l1)
        {
            foreach (var node in l1.Nodes)
            {
                yield return node;
            }
        }
        else if (this is LNCommunity ln)
        {
            foreach (var community in ln.Communities)
            {
                foreach (var node in community.FetchAllNodes())
                {
                    yield return node;
                }
            }
        }
        else
        {
            throw new InvalidOperationException("Unknown community type");
        }
    }
}
public record L1Community(HashSet<int> Nodes) : Community;
public record LNCommunity(Community[] Communities) : Community;
public struct NoEdgeData;

public sealed partial class Graph<N, E>
{
    public Dictionary<int, int> InitialCommunityAssignments()
    {
        var nodeCount = _Nodes.Count;
        var assignments = new Dictionary<int, int>();
        for (int i = 0; i < nodeCount; i++)
        {
            assignments[i] = i;
        }
        return assignments;
    }

    public float ComputeModularity(Dictionary<int, int> communityAssignments, (int NodeId, int Community)? localMove)
    {

        // [formula]
        //  Q = 1/2m * einsum(i,j -> A_ij - k_i * k_j / 2m) * delta(c_i, c_j)
        //  where
        //   A_ij is the weight of the edge between node i and node j
        //   k_i is the degree of node i
        //   delta(c_i, c_j) is 1 if node i and node j are in the same community, otherwise 0
        //   m is the total weight of all edges in the graph
        // [implementation]
        // 1. the coefficient to the final result is not significant during the optimization process
        // 2. use the symmetry of the matrix A_ij to reduce the number of operations

        var m = _AllWeight;
        var nodeCount = _Nodes.Count;
        float score = 0.0f;

        if (localMove is null)
        {
            for (int i = 0; i < nodeCount; i++)
            {
                var assigni = communityAssignments[i];
                for (int j = i + 1; j < nodeCount; j++)
                {
                    var assignj = communityAssignments[j];
                    if (assigni != assignj)
                    {
                        continue;
                    }
                    var ki = _Connections[i].Count;
                    var kj = _Connections[j].Count;
                    if (_Connections[i].TryGetValue(j, out var edgeij))
                    {
                        var edgeIjWeight = _Edges[edgeij].Weight;
                        var scoreDelta = edgeIjWeight - (ki * kj / (m + m));
                        score += scoreDelta;
                    }
                    else
                    {
                        score -= ki * kj / (m + m);
                    }

                }
            }
        }
        else
        {
            var (movedNodeId, movedCommunityTarget) = localMove.Value;
            for (int i = 0; i < nodeCount; i++)
            {
                var assigni = i == movedNodeId ? movedCommunityTarget : communityAssignments[i];
                for (int j = i + 1; j < nodeCount; j++)
                {
                    var assignj = j == movedNodeId ? movedCommunityTarget : communityAssignments[j];
                    if (assigni != assignj)
                    {
                        continue;
                    }
                    var ki = _Connections[i].Count;
                    var kj = _Connections[j].Count;
                    if (_Connections[i].TryGetValue(j, out var edgeij))
                    {
                        var edgeIjWeight = _Edges[edgeij].Weight;
                        var scoreDelta = edgeIjWeight - (ki * kj / (m + m));
                        score += scoreDelta;
                    }
                    else
                    {
                        score -= ki * kj / (m + m);
                    }

                }
            }
        }


        return score / m;
    }

    record struct ProposedMove(int NodeId, int Community, float Modularity);

    record class ParallelLocalData(
        ConcurrentBag<ProposedMove> Moves,
        Graph<N, E> Graph,
        float CurrentModularity,
        Dictionary<int, int> CommunityAssignments
    );


#if DEBUG
    const int ParallelThreshold = 0;
#else
    const int ParallelThreshold = 512;
#endif

    public void OptimizeModularity(Dictionary<int, int> communityAssignments)
    {
        var currentModularity = ComputeModularity(communityAssignments, null);

        float previousModularity;
        var nodeCount = _Nodes.Count;

        if (nodeCount > ParallelThreshold)
        {
#if DEBUG
            Console.WriteLine($"DEBUG: using parallel fast local move");
#endif
            // parallel & batch fast local move
            ConcurrentBag<ProposedMove> batchMoving = [];
            do
            {
                previousModularity = currentModularity;

                static void Finally(ParallelLocalData _) { };
                static ParallelLocalData LoopBody(int i, ParallelLoopState _, ParallelLocalData localData)
                {
                    var node = new NodeId(i);
                    var modularity = localData.CurrentModularity;
                    var newAssign = localData.Graph.FastLocalMove(node, localData.CommunityAssignments, ref modularity);
                    if (newAssign is not null)
                    {
                        localData.Moves.Add(new ProposedMove(i, newAssign.Value, modularity));
                    }
                    return localData;
                };

                var local = new ParallelLocalData(batchMoving, this, currentModularity, communityAssignments);
                Parallel.For(0, nodeCount, () => local, LoopBody, Finally);

                // perform the moves
                // XXX: sort by modularity and use only the top N?
                foreach (var move in batchMoving)
                {
                    communityAssignments[move.NodeId] = move.Community;
                }

                currentModularity = ComputeModularity(communityAssignments, null);

                // clear the bag
                batchMoving.Clear();

            } while (currentModularity > previousModularity);

        }
        else
        {
            // batch fast local move
            List<ProposedMove> batchMoving = [];
            do
            {
                previousModularity = currentModularity;

                for (int i = 0; i < nodeCount; i++)
                {
                    var node = new NodeId(i);
                    var modularity = currentModularity;
                    var newAssign = FastLocalMove(node, communityAssignments, ref modularity);
                    if (newAssign is not null)
                    {
                        batchMoving.Add(new ProposedMove(i, newAssign.Value, modularity));
                    }
                }

                // perform the moves
                // XXX: sort by modularity and use only the top N?
                foreach (var move in batchMoving)
                {
                    communityAssignments[move.NodeId] = move.Community;
                }

                currentModularity = ComputeModularity(communityAssignments, null);

                // clear the moves
                batchMoving.Clear();

            } while (currentModularity > previousModularity);
        }
    }

    int? FastLocalMove(NodeId node, Dictionary<int, int> communityAssignments, ref float currentModularity)
    {
        var neighbors = _Connections[node.Id];
        int bestAssign = communityAssignments[node.Id];
        bool changed = false;
        foreach (var neighbor in neighbors.Keys)
        {
            var neighborAssign = communityAssignments[neighbor];
            if (neighborAssign == bestAssign)
            {
                continue;
            }
            var newModularity = ComputeModularity(communityAssignments, (node.Id, neighborAssign));
            if (newModularity > currentModularity)
            {
                // commit the assignments
                bestAssign = neighborAssign;
                currentModularity = newModularity;
                changed = true;
            }
        }
        if (changed)
        {
            return bestAssign;
        }
        else
        {
            return null;
        }

    }

    public List<HashSet<int>> Refine(Dictionary<int, int> communityAssignments)
    {
        // this is the community assignments by louvain
        // each community might get split into multiple communities
        // if there are partitions that are not connected to each other
        List<HashSet<int>> communitiesByLouvain = [];


        // fill `communitiesByLouvain`
        {
            // maps louvain community to the new community
            Dictionary<int, int> relabel = [];

            int getRelabeledCommunity(int louvainCommunity)
            {
                if (!relabel.TryGetValue(louvainCommunity, out var relabeledCommunity))
                {
                    relabeledCommunity = relabel.Count;
                    Debug.Assert(relabeledCommunity == communitiesByLouvain.Count);
                    relabel[louvainCommunity] = relabeledCommunity;
                    communitiesByLouvain.Add([]);
                }
                return relabeledCommunity;
            }

            foreach (var kv in communityAssignments)
            {
                var louvainCommunity = kv.Value;
                var relabeledCommunity = getRelabeledCommunity(louvainCommunity);
                communitiesByLouvain[relabeledCommunity].Add(kv.Key);
            }
        }

        // XXX: parallelize?
        // validate the inner connections in each community
        for (int i = 0; i < communitiesByLouvain.Count; i++)
        {
            var community = communitiesByLouvain[i];
            if (community.Count == 1)
            {
                continue;
            }

            // copy the community
            var leftMembers = new HashSet<int>(community);
            var queue = new Queue<int>();
            queue.Enqueue(community.First());
            while (queue.TryDequeue(out var node))
            {
                var newlyVisited = leftMembers.Remove(node);
                if (!newlyVisited)
                {
                    // already visited, skip
                    continue;
                }

                var neighbors = _Connections[node];
                foreach (var neighbor in neighbors.Keys)
                {
                    if (!community.Contains(neighbor))
                    {
                        // the subgraph cannot be connected via this node
                        continue;
                    }
                    queue.Enqueue(neighbor);
                }
            }

            if (leftMembers.Count > 0)
            {
                // the subgraph is not connected
                // we need to split the community
                // P.S:
                // due to the nature of this iterative algorithm,
                // we can simply add the split task into the queue

                community.ExceptWith(leftMembers);
                communitiesByLouvain.Add(leftMembers);
            }
        }

        return communitiesByLouvain;
    }


    /// <summary>
    ///  build the hierarchy of communities using leiden algorithm
    /// </summary>
    public List<Graph<Community, NoEdgeData>> Leiden(int? maxIter = null)
    {
        List<Graph<Community, NoEdgeData>> result = [];
        Graph<Community, NoEdgeData> highLevelGraph;
        {
            var g = this;
            var g1 = g.LeidenL1();
            result.Add(g1);
            var g2 = g1.LeidenLN();
            if (g1.CountNodes() == g2.CountNodes())
            {
                return result;
            }
            highLevelGraph = g2;
        }

        var count = highLevelGraph.CountNodes();
        int previous;

        if (maxIter is null)
        {
            int iter = 0;
            do
            {
                previous = count;
                result.Add(highLevelGraph);
                highLevelGraph = highLevelGraph.LeidenLN();
                count = highLevelGraph.CountNodes();
                if (++iter % 100 == 0)
                {
                    Console.WriteLine($"Leiden iteration {iter}: {count}");
                }
            } while (previous != count);
        }
        else
        {

            for (int i = 0; i < maxIter; i++)
            {
                previous = count;
                result.Add(highLevelGraph);
                highLevelGraph = highLevelGraph.LeidenLN();
                count = highLevelGraph.CountNodes();
                if (i % 100 == 0)
                {
                    Console.WriteLine($"Leiden iteration {i}: {count}");
                }
                if (previous == count)
                {
                    break;
                }
            }
        }

        return result;
    }
}

public static class CommunityExtensions
{
    internal static Graph<Community, NoEdgeData> LeidenL1<N, E>(this Graph<N, E> g)
    {
        var initialAssignments = g.InitialCommunityAssignments();
        g.OptimizeModularity(initialAssignments);
        var communities = g.Refine(initialAssignments);
        return g.CompressL1(communities);
    }

    internal static Graph<Community, NoEdgeData> LeidenLN(this Graph<Community, NoEdgeData> g)
    {
        var initialAssignments = g.InitialCommunityAssignments();
        g.OptimizeModularity(initialAssignments);
        var communities = g.Refine(initialAssignments);
        return g.CompressLN(communities);
    }

    internal static Graph<Community, NoEdgeData> CompressL1<N, E>(this Graph<N, E> g, List<HashSet<int>> communities)
    {
        Dictionary<int, int> nodeToCommunity = [];
        var newGraph = new Graph<Community, NoEdgeData>();
        for (int i = 0; i < communities.Count; i++)
        {
            var community = communities[i];
            var newCommunity = new L1Community(community);
            newGraph.AddNode(newCommunity);
            foreach (var node in community)
            {
                nodeToCommunity[node] = i;
            }
        }


        var nodeCount = g.CountNodes();
        for (int i = 0; i < nodeCount; i++)
        {
            var assigni = nodeToCommunity[i];
            for (int j = i + 1; j < nodeCount; j++)
            {
                var assignj = nodeToCommunity[j];
                if (assigni != assignj)
                {
                    continue;
                }

                var edgeInfo = g.GetEdgeInfo(new NodeId(i), new NodeId(j));
                if (edgeInfo is not null)
                {
                    newGraph.AddEdge(
                        new NodeId(assigni),
                        new NodeId(assignj),
                        new NoEdgeData(),
                        edgeInfo.Value.Weight
                    );
                }
            }
        }

        return newGraph;
    }

    internal static Graph<Community, NoEdgeData> CompressLN(this Graph<Community, NoEdgeData> g, List<HashSet<int>> communities)
    {
        Dictionary<int, int> nodeToCommunity = [];
        var newGraph = new Graph<Community, NoEdgeData>();
        for (int i = 0; i < communities.Count; i++)
        {
            var community = communities[i];
            var newCommunity = new LNCommunity(community.Select(c => g.GetNodeData(new NodeId(c))).ToArray());
            newGraph.AddNode(newCommunity);
            foreach (var node in community)
            {
                nodeToCommunity[node] = i;
            }
        }


        var nodeCount = g.CountNodes();
        for (int i = 0; i < nodeCount; i++)
        {
            var assigni = nodeToCommunity[i];
            for (int j = i + 1; j < nodeCount; j++)
            {
                var assignj = nodeToCommunity[j];
                if (assigni != assignj)
                {
                    continue;
                }

                var edgeInfo = g.GetEdgeInfo(new NodeId(i), new NodeId(j));
                if (edgeInfo is not null)
                {
                    newGraph.AddEdge(
                        new NodeId(assigni),
                        new NodeId(assignj),
                        new NoEdgeData(),
                        edgeInfo.Value.Weight
                    );
                }
            }
        }

        return newGraph;
    }
}

public static partial class Test
{
    public static void TestAlgo()
    {
        List<(string, string, float)> edges = [
            // from | to | weight
            ("Fortran", "C", 0.5f),
            ("Fortran", "LISP", 0.3f),
            ("Fortran", "MATLAB", 0.6f),
            ("C", "C++", 0.9f),
            ("C", "Java", 0.2f),
            ("C", "Go", 0.6f),
            ("LISP", "ML", 0.5f),
            ("LISP", "OCaml", 0.2f),
            ("LISP", "Haskell", 0.2f),
            ("LISP", "Ruby", 0.5f),
            ("LISP", "Julia", 0.6f),
            ("ML", "OCaml", 0.8f),
            ("ML", "Haskell", 0.5f),
            ("OCaml", "Haskell", 0.3f),
            ("OCaml", "F#", 0.6f),
            ("Haskell", "Julia", 0.2f),
            ("C++", "Java", 0.5f),
            ("C++", "Python", 0.32f),
            ("C++", "Ruby", 0.2f),
            ("C++", "C#", 0.5f),
            ("Java", "Ruby", 0.4f),
            ("Java", "Python", 0.5f),
            ("Java", "C#", 0.6f),
            ("Java", "Go", 0.45f),
            ("Java", "Julia", 0.1f),
            ("Python", "F#", 0.2f),
            ("Python", "Julia", 0.4f),
            ("C#", "F#", 0.3f),
        ];
        Dictionary<string, NodeId> nodes = new();
        var g = new Graph<string, string>();
        foreach (var (from, to, weight) in edges)
        {
            if (!nodes.TryGetValue(from, out var fromId))
            {
                fromId = nodes[from] = g.AddNode(from);
            }
            if (!nodes.TryGetValue(to, out var toId))
            {
                toId = nodes[to] = g.AddNode(to);
            }
            g.AddEdge(fromId, toId, $"{from}->{to}", weight);
        }

        var hierarchy = g.Leiden();
        var gLast = hierarchy.Last();
        for (int i = 0; i < gLast.CountNodes(); i++)
        {
            Console.WriteLine($"Community {i}:");
            foreach (var node in gLast.GetNodeData(new NodeId(i)).FetchAllNodes())
            {
                var nodInOriinalGraph = new NodeId(node);
                Console.WriteLine($"  {g.GetNodeData(nodInOriinalGraph)}");
            }
        }
    }
}