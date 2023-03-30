
//
// Created by Tianyu Zhang on 3/15/23.
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//

#include "klc3/FlowAnalysis/FlowGraphVisualizer.h"
#include "klc3/FlowAnalysis/BackBoneAnalyzer.h"
#include "klc3/Core/MemoryValue.h"
#include "klc3/Core/State.h"

namespace klc3 {

Node *SCCSet::getSCCRoot(Node *node) const {
    auto it = info.find(node);
    assert(it != info.end() && "SCC root not determined");
    // Can also be outside the current subroutine
    // e.g. RET_EDGE to the caller
    // if (it != info.end()) {
    //     return nullptr;
    // }

    return it->second.sccRoot;
}

void SCCSet::analysisSCC() {
    dfsClock = 0;
    assert(tarjanStack.empty() && "Tarjan stack is not empty");
    sccBlockRoot.clear();
    info.clear();

    // Backtrack DFS
    for (auto &node : nodes) {
        if (info.find(node) == info.end()) {  // not visited
            tarjanDFS(node);
        }
    }
}

void SCCSet::tarjanDFS(Node *v) {
    ++dfsClock;

    info[v] = {
            dfsClock,
            dfsClock,
            nullptr
    };  // also marked as visited
    SCCSet::dfsInfo &vInfo = info[v];

    tarjanStack.push(v);

    Node *w;
    for (auto &e : v->allOutEdges()) {
        if (e->type() == Edge::JSR_EDGE || e->type() == Edge::JSRR_EDGE) {
            continue;  // do not analyze across into subroutines, but recognize RET
        }
        // Use SUBROUTINE_VIRTUAL_EDGE in place of JSR_EDGE (there is really no JSRR_EDGE during static analysis stage).

        // Ignore node virtual edges
        if (e->type() == Edge::LOOP_H2H_EDGE || e->type() == Edge::LOOP_H2X_EDGE) continue;

        // Ignore likely uncoverable edges
        if (e->flags & Edge::LIKELY_UNCOVERABLE) continue;

        w = e->to();
        if (w == nullptr || !containNode(w)) continue;

        if (info.find(w) == info.end()) {
            tarjanDFS(w);
            vInfo.dfsLow = std::min(vInfo.dfsLow, info[w].dfsLow);
        } else if (info[w].sccRoot == nullptr) {
            vInfo.dfsLow = std::min(vInfo.dfsLow, info[w].dfsPre);
        }
    }

    if (vInfo.dfsLow == vInfo.dfsPre) {
        int sccNodeCount = 0;
        do {
            w = tarjanStack.top();
            tarjanStack.pop();
            info[w].sccRoot = v;
            ++sccNodeCount;
        } while (w != v);

        if (sccNodeCount > 1) {
            sccBlockRoot.emplace(v);
        } else {
            sccSingle.emplace(v);
        }

        sccRoots.emplace(v);
    }
}

// bool FlowGraphVisualizer::shouldStopContinuousBlock(VisualNode *node) const {
//     assert(node->outEdges.size() == 1 && "Not a continuous block");
//     if (node->notMergeAfterward) return true;                     // current node requires not to merge
//     VisualNode *next = node->outEdges[0]->to;
//     if (next->inEdges.size() >= 2) return true;                   // next node has multiple entries
//     if (!next->label.empty()) return true;                        // next node has label
//     if (next->blockColors != node->blockColors) return true;      // next node has different color
//     if (next->contentColor != node->contentColor) return true;    // next node has different color
//     if (next->notMergeAhead) return true;                         // next node requires not to merge
//     if (node->outEdges[0]->notMerge) return true;                 // unique outEdge requires not to merge
//     if (node->visible != next->visible) return true;              // different visibility
//     if (nodeParent.find(node)->second != nodeParent.find(next)->second) return true;  // different layer
//     return false;
// }


bool BackBoneAnalyzer::eliminateCycles(const SCCSet &ss) {

}

bool BackBoneAnalyzer::compressBackBone(Node *curNode, BackBoneNode *backBoneNode) {

    llvm::outs() << curNode->name() << "\n";

    // Add current node into backbone
    // backBoneNode->sccRoots_.emplace(curNode);
    // backBoneNode->nodes_.emplace(curNode);
    backBoneNode->sg.nodes.emplace(curNode);
    visited.emplace(curNode, backBoneNode);

    vector<Edge *> outEdges = curNode->allOutEdges();
    // InstValue::InstID id = curNode->inst().get()->instID();

    for (auto &edge : outEdges) {
        llvm::outs() << "\t" << "edge" << "\n";
        if (edge->from() != nullptr)
            llvm::outs() << "\t" << edge->from()->name() << "\n";
        if (edge->to() != nullptr)
            llvm::outs() << "\t" << edge->to()->name() << "\n";
        llvm::outs() << "\t" << edge->type() << "\n";

        // 
        // For now, we ignore subroutines as we have not implemented subroutine expansion yet
        // TODO: add subroutine expansion here
        // 
        if (edge->type() == Edge::JSR_EDGE || edge->type() == Edge::JSRR_EDGE || edge->type() == Edge::RET_EDGE)
            continue;  // do not analyze across into subroutines

        // We will go for SUBROUTINE_VIRTUAL_EDGE, assuming the subroutine returns properly
        // Just treat this as a NORMAL_EDGE

        if (edge->flags & Edge::LIKELY_UNCOVERABLE) continue;

        // Ignore node virtual edges
        if (edge->type() == Edge::LOOP_H2H_EDGE || edge->type() == Edge::LOOP_H2X_EDGE) continue;

        // Last node to exit the subroutine or program
        if (edge->to() == nullptr) continue;

        // 
        // Right now, if we see nodes already examined,
        // we stop and treat the dest_node as in a different backbone node
        // 
        Node *dest_node = edge->to();
        if (visited.find(dest_node) != visited.end()) {
            if (visited.find(dest_node) != visited.find(curNode))
                backBoneNode->outEdges_.emplace(edge);
            continue;
        }

        // Check if the dest_node has multiple incoming edges
        if (dest_node->allInEdges().size() > 1) {
            backBoneNode->outEdges_.emplace(edge);
            analyzeBackBone(dest_node);
            continue;
        }
s        
        // Check if we meet BR
        // This is where we make new branches &
        // Check consecutive BRs
        
        if (edge->flags & Edge::Flags::BR_BRANCH_EDGE) {
            backBoneNode->condition_flag_ = BackBoneNode::ConditionFlags(curNode->inst().get()->cc() | backBoneNode->condition_flag_);

            if (backBoneNode->condition_flag_ == BackBoneNode::ConditionFlags::BRnzp) {
                compressBackBone(dest_node, backBoneNode);
                continue;
            }

            backBoneNode->outEdges_.emplace(edge);
            analyzeBackBone(dest_node);
            continue;
        }

        if (edge->flags & Edge::Flags::BR_CONTINUE_EDGE) {
            backBoneNode->condition_flag_ = BackBoneNode::ConditionFlags(curNode->inst().get()->cc() | backBoneNode->condition_flag_);

            if (backBoneNode->condition_flag_ == BackBoneNode::ConditionFlags::BRnzp)
                continue;

            // Check next instruction is BR
            if (dest_node->inst().get()->instID() == InstValue::BR) {
                compressBackBone(dest_node, backBoneNode);
                continue;
            }

            backBoneNode->outEdges_.emplace(edge);
            analyzeBackBone(dest_node);
            continue;
        }

        compressBackBone(dest_node, backBoneNode);
    }
    return true;
}


// 
// TODO: talk about multiple incoming and outcoming edges !!!!!
// 
// bool BackBoneAnalyzer::analyzeBackBoneSCCSingleDFS(Node *curNode, BackBoneNode *backBoneNode, const SCCSet &ss) {
//     llvm::outs() << curNode->name() << "\n";

//     // Add current node into backbone
//     backBoneNode->sccRoots_.emplace(curNode);
//     // backBoneNode->nodes_.emplace(curNode);
//     backBoneNode->sg.nodes.emplace(curNode);
//     visited.emplace(curNode, backBoneNode);

//     vector<Edge *> outEdges = curNode->allOutEdges();
//     InstValue::InstID id = curNode->inst().get()->instID();
//     for (auto &edge : outEdges) {

//         if (edge->from() != nullptr)
//             llvm::outs() << "\t" << edge->from()->name() << "\n";
//         if (edge->to() != nullptr)
//             llvm::outs() << "\t" << edge->to()->name() << "\n";
//         llvm::outs() << "\t" << edge->type() << "\n";

//         // 
//         // For now, we ignore subroutines as we have not implemented subroutine expansion yet
//         // TODO: add subroutine expansion here
//         // 
//         if (edge->type() == Edge::JSR_EDGE || edge->type() == Edge::JSRR_EDGE || edge->type() == Edge::RET_EDGE)
//             continue;  // do not analyze across into subroutines

//         // We will go for SUBROUTINE_VIRTUAL_EDGE, assuming the subroutine returns properly
//         // Just treat this as a NORMAL_EDGE

//         if (edge->flags & Edge::LIKELY_UNCOVERABLE) continue;

//         // Ignore node virtual edges
//         if (edge->type() == Edge::LOOP_H2H_EDGE || edge->type() == Edge::LOOP_H2X_EDGE) continue;

//         // Last node to exit the subroutine or program
//         if (edge->to() == nullptr) continue;

//         // 
//         // Right now, if we see nodes already examined,
//         // we stop and treat the dest_node as in a different backbone node
//         // 
//         Node *dest_node = edge->to();
//         if (visited.find(dest_node) != visited.end()) {
//             backBoneNode->outEdges_.emplace(edge);
//             continue;
//         }

//         // Check if the dest_node has multiple incoming edges
//         if (dest_node->allInEdges().size() > 1) {
//             // analyzeBackBone(dest_node, ss);
//             continue;
//         }

//         // 
//         // Check if we meet BR
//         // This is where we make new branches
//         // TODO:
//         // 1. Check consecutive BRs
//         // 2. Eliminate "BR LABEL" and "BRnzp LABEL"
//         // !!!!
//         // 3. DFS or BFS???
//         // 
//         if (id == InstValue::BR) {
//             backBoneNode->outEdges_.emplace(edge);
//             // analyzeBackBone(dest_node, ss);
//             continue;
//         }

//         // 
//         // Check if we move into a scc with multiple nodes
//         // 
//         Node *nextSCCRoot = ss.getSCCRoot(dest_node);
//         if (ss.sccBlockRoot.find(nextSCCRoot) != ss.sccBlockRoot.end()) {
//             analyzeBackBoneSCCMultipleDFS(dest_node, nextSCCRoot, backBoneNode, ss);
//         } else if (ss.sccSingle.find(nextSCCRoot) != ss.sccSingle.end()) {
//             analyzeBackBoneSCCSingleDFS(dest_node, backBoneNode, ss);
//         } else {
//             // Should not reach here
//         }
//     }

//     return true;
// }


// bool BackBoneAnalyzer::analyzeBackBoneSCCMultipleDFS(Node *curNode, Node *sccRootNode, BackBoneNode *backBoneNode, const SCCSet &ss) {

//     llvm::outs() << curNode->name() << "\n";

//     // Check if we have gone through a cycle
//     // if (backBoneNode->nodes().find(curNode) != backBoneNode->nodes().end())
//     if (backBoneNode->subgraph().containNode(curNode))
//         return true;
//     // Add current node into backbone
//     if (backBoneNode->sccRoots().find(sccRootNode) == backBoneNode->sccRoots().end())
//         backBoneNode->sccRoots_.emplace(sccRootNode);
//     // backBoneNode->nodes_.emplace(curNode);
//     backBoneNode->sg.nodes.emplace(curNode);
//     visited.emplace(curNode, backBoneNode);

//     // 
//     // We do not care about the instrcution type in a scc with multiple nodes
//     // We only care when we get nodes outside of scc
//     // TODO: 
//     // 1. MERGE node if we see only one exit edge in this scc to the outside
//     // 2. Need to decide what to do if we see more than one exit edges
//     // 
//     // InstValue::InstID id = curNode->inst().get()->instID();
//     // 
    
//     // Go through all the edges of current node
//     vector<Edge *> outEdges = curNode->allOutEdges();
//     llvm::outs() << "\tnumber:" << outEdges.size() << "\n";

//     bool has_go_out = false;
//     for (auto &edge : outEdges) {

//         if (edge->from() != nullptr)
//             llvm::outs() << "\t" << edge->from()->name() << "\n";
//         if (edge->to() != nullptr)
//             llvm::outs() << "\t" << edge->to()->name() << "\n";
//         llvm::outs() << "\t" << edge->type() << "\n";

//         // 
//         // For now, we ignore subroutines as we have not implemented subroutine expansion yet
//         // TODO: add subroutine expansion here
//         // 
//         if (edge->type() == Edge::JSR_EDGE || edge->type() == Edge::JSRR_EDGE || edge->type() == Edge::RET_EDGE)
//             continue;  // do not analyze across into subroutines

//         // We will go for SUBROUTINE_VIRTUAL_EDGE, assuming the subroutine returns properly
//         // Just treat this as a NORMAL_EDGE

//         if (edge->flags & Edge::LIKELY_UNCOVERABLE) continue;

//         // Ignore node virtual edges
//         if (edge->type() == Edge::LOOP_H2H_EDGE || edge->type() == Edge::LOOP_H2X_EDGE) continue;

//         // Last node to exit the subroutine or program
//         if (edge->to() == nullptr) continue;

//         // Check if the edge goes outside the scc
//         Node *dest_node = edge->to();
//         Node *newRootNote = ss.getSCCRoot(dest_node);

//         // Check if the dest_node has multiple incoming edges
//         if (dest_node->allInEdges().size() > 1) {
//             // analyzeBackBone(dest_node, ss);
//             continue;
//         }

//         if (newRootNote == sccRootNode) {  // still inside the scc
//             analyzeBackBoneSCCMultipleDFS(dest_node, sccRootNode, backBoneNode, ss);
//             continue;
//         }
        
//         // 
//         // goes outside the current scc
//         // We want to merge nodes if this is the first time to go outside current scc
//         // If not, we treat the dest_node as in a different backbone node
//         // 
//         if (visited.find(dest_node) != visited.end()) {
//             backBoneNode->outEdges_.emplace(edge);
//             continue;
//         }

//         // analyzeBackBone(dest_node, ss);

//         // if (has_go_out) {
//         //     analyzeBackBone(dest_node, ss);
//         //     continue;
//         // }
//         // // 
//         // // Check if we move into a scc with multiple nodes
//         // // 
//         // has_go_out = true;
//         // Node *nextSCCRoot = ss.getSCCRoot(dest_node);
//         // if (ss.sccBlockRoot.find(nextSCCRoot) != ss.sccBlockRoot.end()) {
//         //     analyzeBackBoneSCCMultipleDFS(dest_node, nextSCCRoot, backBoneNode, ss);
//         // } else if (ss.sccSingle.find(nextSCCRoot) != ss.sccSingle.end()) {
//         //     analyzeBackBoneSCCSingleDFS(dest_node, backBoneNode, ss);
//         // } else {
//         //     // Should not reach here
//         // }
//     }

//     return true;
// }


// 
// When first called,
// entryNode must be the entryNode of the subroutine
// Returns false if meet spaghetti code
// 
bool BackBoneAnalyzer::analyzeBackBone(Node *entryNode) {
    llvm::outs() << "new\n";

    // 
    // One BackBoneNode can have multiple sccs, and all BackBoneNodes form a 
    // singlely-connected graph (Tree), so we need to analyze from the root, which is the entry node
    // 
    // In fact, a single backbone node contains the set of sccs starting from the root of a subtree,
    // ending with the place where we have a complete branch (flags = nzp).
    // e.g.
    /*
        NODE1
        NODE2   --->    NODE5
        NODE3
        NODE4

        Each node above is a scc (no loops here)
        We have three backbone nodes:
        BackBoneNode1: NODE1, NODE2
        BackBoneNode2: NODE3, NODE4
        BackBoneNode3: NODE5

        Tree structure shown below:

              1
            /   \
           2     3
    */

    // Entry must not in the backbone list yet
    // assert(backBoneNodes.find(entryNode) == backBoneNodes.end() && "Duplicated BackBoneNode");
    if (backBoneNodes.find(entryNode) != backBoneNodes.end()) {
        return true;
    }
    

    // Manage all nodes in a tree based structure starting from the root

    // 
    // Construct new backBoneNode information structure
    // 
    BackBoneNode *backBoneNode = new BackBoneNode;
    backBoneNode->name_ = entryNode->inst()->getLabel();
    backBoneNode->entryNode_ = entryNode;
    // backBoneNode->sccRoots_.emplace(entryNode);
    backBoneNodes.emplace(entryNode, backBoneNode);
    backBoneNode->la = new LoopAnalyzer(fg);            // Might be buggy, not sure...............
    visited.emplace(entryNode, backBoneNode);

    // If the root for backbone
    if (entryNode == backBoneEntryPoint.first)
        backBoneEntryPoint.second = backBoneNode;


    compressBackBone(entryNode, backBoneNode);


    // 
    // Add all nodes in the current scc into the backboneNode
    // Then move along the path
    // 
    // Node *currentSCCRoot = ss.getSCCRoot(entryNode);
    // if (ss.sccBlockRoot.find(currentSCCRoot) != ss.sccBlockRoot.end()) {
    //     // backBoneNode->nodes_.emplace(node);
    //     analyzeBackBoneSCCMultipleDFS(currentSCCRoot, currentSCCRoot, backBoneNode, ss);
    // } else if (ss.sccSingle.find(currentSCCRoot) != ss.sccSingle.end()) {
    //     // backBoneNode->nodes_.emplace(currentSCCRoot);
    //     analyzeBackBoneSCCSingleDFS(currentSCCRoot, backBoneNode, ss);

    // } else {
    //     // Should not reach here
    // }

    return true;
}


bool BackBoneAnalyzer::analyzeBackBoneNodesLoop() {
    
    bool ret = true;
    for (auto &node : getAllBackBoneNodes()) {
        ret &= node->la->analyzeLoops(node->entryNode(), node->subgraph());
    }
    return ret;
}


set<BackBoneNode *> BackBoneAnalyzer::getAllBackBoneNodes() const {
    set<BackBoneNode *> ret;
    for (auto it : backBoneNodes) {
        ret.emplace(it.second);
    }
    return ret;
}

BackBoneAnalyzer::~BackBoneAnalyzer() {
    for (auto &it : backBoneNodes) delete it.second;
}

void BackBoneAnalyzer::dump(llvm::raw_ostream &os) const {
    int h2xCount = 0;

    dumpBackBoneNodes(os, "");

    os << "Total node count: " << backBoneNodes.size() << "\n";

    for (const auto &it : backBoneNodes) {
        h2xCount += it.second->OutEdges().size();
    }
    os << "Total H2X count: " << h2xCount << "\n";
}

void BackBoneAnalyzer::dumpBackBoneNodes(llvm::raw_ostream &os, const string &indent) const {

    for (const auto &node : backBoneNodes) {

        os << indent << "-> " << node.second->name() << "\n";
        os << indent << "    Entry point: " << node.second->entryNode()->addr() << "\n";
        os << indent << "    Entry node name: " << node.second->entryNode()->name() << "\n";
        os << indent << "    Num of nodes count: " << node.second->subgraph().nodes.size() << "\n";
        os << indent << "    Num of scc Nodes count: " << node.second->sccRoots().size() << "\n";
        // os << indent << "    H2X segment count: " << node.second->h2xSegments().size() << "\n";
        os << indent << "    H2X edge count: " << node.second->OutEdges().size() << "\n";
    }
}

}
