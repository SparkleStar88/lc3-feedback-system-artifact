
//
// Created by Tianyu Zhang on 3/15/23.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//

#ifndef KLC3_BACKBONEANALYZER_H
#define KLC3_BACKBONEANALYZER_H

#include "FlowGraph.h"
#include "LoopAnalyzer.h"

class Node;
class Edge;

namespace klc3 {

class SCCSet : public Subgraph {
public:

    // explicit SCCSet(Subgraph sg) : nodes(sg.nodes) {}
    explicit SCCSet(Subgraph sg) { nodes = sg.nodes; }

    struct dfsInfo {
        int dfsPre;
        int dfsLow;
        Node *sccRoot = nullptr;
    };

    Node *getSCCRoot(Node *node) const;

    void analysisSCC();

    void tarjanDFS(Node *v);

protected:

    friend class BackBoneNode;
    friend class BackBoneAnalyzer;


    unordered_map<Node *, dfsInfo> info;

    unordered_set<Node *> sccBlockRoot;  // root of scc that contains more than one nodes

    unordered_set<Node *> sccSingle;    // node of scc that contains only one node

    unordered_set<Node *> sccRoots;    // root of sccs

    int dfsClock = 0;

    stack<Node *> tarjanStack;
};


class BackBoneNode {
public:

    const string &name() const { return name_; }  // use entrypoint first label or pseudo label

    Node *entryNode() const { return entryNode_; }

    // const auto &nodes() const { return nodes_; }

    const auto &subgraph() const { return sg; }

    const auto &sccRoots() const { return sccRoots_; }

    const set<Path> &Paths() const { return paths_; }

    const set<Edge *> &OutEdges() const { return outEdges_; }

    // const auto &segmentLastEdges() const { return segmentLastEdges_; }  // for updating loop path in PruningSearcher

    // TODO:
    // New: we want to do loop analysis on each backbone node
    LoopAnalyzer *la;

protected:
    friend class BackBoneAnalyzer;

    BackBoneNode() = default;  // can only be created by BackBoneAnalyzer

    string name_;  // use entrypoint first label or pseudo label

    Node *entryNode_ = nullptr;

    // Change: we use Subgraph here since we would like to analyze loops in each structure
    Subgraph sg;

    // All nodes in this backbone node
    // The nodes probably form paths from entry to exit
    // unordered_set<Node *> nodes_;

    // root of sccs that contains in this backbone node
    // including multiple and single nodes scc
    unordered_set<Node *> sccRoots_;

    set<Path> paths_;

    // Edges to other BackBoneNodes
    set<Edge *> outEdges_;
    
    enum ConditionFlags {
        BRnone = 0,
        BRn = 4,
        BRz = 2,
        BRp = 1,
        BRnz = BRn | BRz,
        BRzp = BRz | BRp,
        BRnp = BRn | BRp,
        BRnzp = BRn | BRz | BRp,
        BR = BRnzp
    };

    ConditionFlags condition_flag_ = ConditionFlags::BRnone;

    // unordered_set<Edge *> segmentLastEdges_;
};


class BackBoneAnalyzer {
public:

    explicit BackBoneAnalyzer(FlowGraph *fg) : fg(fg) {}

    BackBoneAnalyzer(const BackBoneAnalyzer &) = delete;

    ~BackBoneAnalyzer();

    // Analyze backbone for each subroutine
    // In klc3.cpp, now we do not have subroutine expansion,
    // so we do the analysis subroutine-wise.
    bool analyzeBackBone(Node *entryNode);

    // Used to analyze scc with one node
    // bool analyzeBackBoneSCCSingleDFS(Node *curNode, Node *backBoneEntryNode, BackBoneNode *backBoneNode, const SCCSet &sss);
    bool analyzeBackBoneSCCSingleDFS(Node *curNode, BackBoneNode *backBoneNode, const SCCSet &ss);

    // Used to analyze scc with more than one node
    bool analyzeBackBoneSCCMultipleDFS(Node *curNode, Node *sccRootNode, BackBoneNode *backBoneNode, const SCCSet &ss);

    // 
    // 
    // New
    // 
    bool compressBackBone(Node *curNode, BackBoneNode *backBoneNode);

    bool eliminateCycles(const SCCSet &ss);


    // Analyze loops of all backbone nodes
    bool analyzeBackBoneNodesLoop();

    // ================== Helper functions ==================== //
    const unordered_map<Node *, BackBoneNode *> &getBackBoneVisitedMap() const { return visited; }

    pair<Node*, BackBoneNode*> &getBackBoneEntryPoint() { return backBoneEntryPoint; }

    set<BackBoneNode*> getAllBackBoneNodes() const;

    void dump(llvm::raw_ostream &os) const;

private:

    FlowGraph *fg;

    // ================= Generalized BackBoneNode Detection ======================== //

    // Map the node to the backbone node that it is contained in
    pair<Node*, BackBoneNode*> backBoneEntryPoint;  // Entry of the backbone structure

    unordered_map<Node *, BackBoneNode *> visited;   // All nodes -> BackBoneNode

    unordered_map<Node*, BackBoneNode *> backBoneNodes;  // entry node -> BackBoneNode

    // For constructing segments, do not use now
    stack<Node *> loopHeadStack;

    void dumpBackBoneNodes(llvm::raw_ostream &os, const string &indent) const;
};

}

#endif //KLC3_BACKBONEANALYZER_H
