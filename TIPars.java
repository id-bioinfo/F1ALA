import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.BitSet;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Queue;
import java.util.Set;
import java.util.Map.Entry;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.CopyOnWriteArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ForkJoinPool;
import java.util.concurrent.ForkJoinTask;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicReference;
import java.util.concurrent.locks.ReentrantLock;
import dr.app.tools.NexusExporter;
import dr.evolution.io.Importer.ImportException;
import dr.evolution.io.NewickImporter;
import dr.evolution.tree.FlexibleNode;
import dr.evolution.tree.FlexibleTree;
import dr.evolution.tree.Tree;
import dr.evolution.util.Taxon;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.nio.file.StandardCopyOption;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.stream.Stream;
import java.nio.file.StandardOpenOption;
import java.nio.charset.StandardCharsets;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;



public class TIPars {

	private String OUTPUT_FOLDER;
	private static boolean isMultiplePlacements = true; /// is output multiple placements
	private static boolean isFastaFile = true; /// is fasta file
	
	//parallelStream setting
	private static int THREAD_NUM = 0; //new add on 20230324
	private static ForkJoinPool forkJoinPool = null;

	private boolean DEBUG = false; /// for debug
	private boolean OUTPUT_PNODE = false; /// is output P node sequence

	private static Tree mytree = null;
	private static String internalnode_nidname = "label";
	private static HashMap<FlexibleNode, Integer> node2edge = new HashMap<FlexibleNode, Integer>();
	private static String[] placements = null;
	private static int edge_number = 1;

	private static ArrayList<String> stringSequencesList = new ArrayList<String>(); /// fasta: store taxaseq and ancseq, ordered follow the seqIdxMap
	private static HashMap<String, Integer> seqIdxMap = new HashMap<String, Integer>(); /// sequence names (taxaseq and ancseq) map to their reading orders from file
	private static int sequence_character_length = -1; /// alignment length
	private static HashMap<FlexibleNode, String> node2seqName = new HashMap<FlexibleNode, String>(); /// tree_nodes map to their sequences name
	private static ArrayList<ConcurrentHashMap<Integer, Byte>> multationSequencesMap = new ArrayList<ConcurrentHashMap<Integer, Byte>>(); /// vcf: store taxaseq and ancseq,ordered follow the seqIdxMap
	private static byte[] ref_sequence = new byte[100000]; /// reference sequence string

	private static double minGlobalMemoryBranchScore = Double.MAX_VALUE;
	private ArrayList<FlexibleNode> minGlobalMemoryBranchScoreNodeList = new ArrayList<FlexibleNode>();

	private static HashMap<String, ArrayList<String>> node2placementSeqName = new HashMap<String, ArrayList<String>>(); /// tree_nodes map to their sequences name for graft (subtree) and assembly
	
	public static char[] alphabet_nt = { 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N', '-' }; /// IUPAC nucleotide codes
	public static char[] alphabet_aa =  {'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', 'B', 'Z', 'X', '-'};  ///IUPAC nucleotide codes
	public static HashMap<Byte, HashSet<Byte>> _nucleotide_nomenclature = null;
	private static double[][] _nucleotide_nomenclature_scoreTable = null;
	private static HashMap<BitSet, Byte> _nucleotide_nomenclature_map2char = null;
	private static double[][] _aminoacid_scoreTable = null;
    private static double[][] _used_scoreTable = null;

	private static ReentrantLock lock = new ReentrantLock();

	public static double MinDoubleNumLimit = 10 * Double.MIN_VALUE;

	public TIPars() {
	}

	// initialization
	public TIPars(Tree mytree, String otype, String output_folder) {
		this.mytree = mytree;
		this.OUTPUT_FOLDER = output_folder;
		setupHashtableOfnode2seqName();
		setupHashtableOfNode2edge(otype);
		node2placementSeqName.clear();
	}

	private void setupHashtableOfNode2edge(String otype) {
		node2edge.clear();
		edge_number = 1;
		if (otype.equals("placement")) {
			for (int i = 0; i < mytree.getExternalNodeCount(); i++) {
				Integer edge = new Integer(edge_number);
				FlexibleNode node = (FlexibleNode) mytree.getExternalNode(i);
				node2edge.put(node, edge);
				edge_number++;
			}

			for (int i = 0; i < mytree.getInternalNodeCount(); i++) {
				Integer edge = new Integer(edge_number);
				FlexibleNode node = (FlexibleNode) mytree.getInternalNode(i);
				node2edge.put(node, edge);
				edge_number++;
			}
		}
	}

	private void setupHashtableOfnode2seqName() {
		node2seqName.clear();
		for (int i = 0; i < mytree.getInternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) mytree.getInternalNode(i);
			String sequenceName = (String) (n.getAttribute(this.internalnode_nidname));
			if (seqIdxMap.containsKey(sequenceName) || seqIdxMap.isEmpty()) {
				node2seqName.put(n, sequenceName);
			} else {
				System.out.println("internalnode=" + sequenceName + " not found in alignment.");
			}
		}
		for (int i = 0; i < mytree.getExternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) mytree.getExternalNode(i);
			String sequenceName = n.getTaxon().getId();
			if (seqIdxMap.containsKey(sequenceName) || seqIdxMap.isEmpty()) {
				node2seqName.put(n, sequenceName);
			} else {
				System.out.println("externalnode=" + sequenceName + " not found in alignment.");
			}
		}
	}

	public static HashMap<String, FlexibleNode> setupHashtableOfseqName2node(Tree tree) {
		HashMap<String, FlexibleNode> mySeqName2node = new HashMap<String, FlexibleNode>();
		for (int i = 0; i < tree.getInternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) tree.getInternalNode(i);
			String sequenceName = (String) (n.getAttribute(internalnode_nidname));
			mySeqName2node.put(sequenceName, n);
		}
		for (int i = 0; i < tree.getExternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) tree.getExternalNode(i);
			String sequenceName = n.getTaxon().getId();
			mySeqName2node.put(sequenceName, n);
		}
		return mySeqName2node;
	}
    
	//new add on 20230221
	public static HashSet<String> setupHashSetOfTaxaName(Tree tree) {
		HashSet<String> mySeqNames = new HashSet<String>();
		for (int i = 0; i < tree.getExternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) tree.getExternalNode(i);
			mySeqNames.add(n.getTaxon().getId());
		}
		return mySeqNames;
	}
	
	private static HashMap<FlexibleNode, String> setupNode2seqName(Tree mytree) {
		HashMap<FlexibleNode, String> mynode2seqName = new HashMap<FlexibleNode, String>();
		for (int i = 0; i < mytree.getInternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) mytree.getInternalNode(i);
			String sequenceName = (String) (n.getAttribute(internalnode_nidname));
			mynode2seqName.put(n, sequenceName);
		}
		for (int i = 0; i < mytree.getExternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) mytree.getExternalNode(i);
			String sequenceName = n.getTaxon().getId();
			mynode2seqName.put(n, sequenceName);
		}
		return mynode2seqName;
	}
	
	// Interface API
	public Tree addQuerySequence(String qname, ConcurrentHashMap<Integer, Byte> nodeQseq, String qid, String pid,
			boolean printDisInfoOnScreen, double[] ABQ_brlen, String otype, int ii) {
		return addQuerySequence_global(qname, nodeQseq, qid, pid, printDisInfoOnScreen, ABQ_brlen, otype, ii);
	}

	public Tree addQuerySequence(String qname, String nodeQseq, String qid, String pid, boolean printDisInfoOnScreen,
			double[] ABQ_brlen, String otype, int ii) {
		return addQuerySequence_global(qname, nodeQseq, qid, pid, printDisInfoOnScreen, ABQ_brlen, otype, ii);
	}

	// Travel all internal nodes for all possible nodeA-nodeB pairs
	// for vcf file
	public Tree addQuerySequence_global(String qname, ConcurrentHashMap<Integer, Byte> nodeQseq, String qid, String pid,
			boolean printDisInfoOnScreen, double[] ABQ_brlen, String otype, int ii) {

		minGlobalMemoryBranchScore = Double.MAX_VALUE;
		minGlobalMemoryBranchScoreNodeList.clear();
		/// if parallel computes branchscore is bigger than minGlobalMemoryBranchScore,
		/// the calculation will stop.
		ArrayList<Integer> nodeIdxAndScoreList = new ArrayList<Integer>(mytree.getNodeCount());
		//// initialize the nodeIdxAndScoreList to store nodeIdx
		for (int i = 0; i < mytree.getNodeCount(); ++i) {
			nodeIdxAndScoreList.add(i);
		}

		//// parallel to compute branchScore for every branch
		ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
			nodeIdxAndScoreList.parallelStream().forEach(nodeIdx -> {
				double score = computeBranchScore(nodeIdx, nodeQseq, qname);
			});
		});
		try {
			forkJoinTask.get();
		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ArrayList<FlexibleNode> selectedNodeList = minGlobalMemoryBranchScoreNodeList;

		//// try to remove ambiguous results
		selectedNodeList = reduceAmbiguousBranch(selectedNodeList, isMultiplePlacements);

		if (DEBUG)
			System.out.println(
					"minQScore: " + minGlobalMemoryBranchScore + ", selectedNodeList:" + selectedNodeList.size());

		MyFlexibleTree best_mynewTree = null;
		String placementStrings = "";
		/// the first element of selectedNodeList is the best node to be inserted
		for (int k = 0; k < selectedNodeList.size(); ++k) {
			/// get sequences
			FlexibleNode selectedNode = selectedNodeList.get(k);
			int selectedNodeBIndex = selectedNode.getNumber();
			ConcurrentHashMap<Integer, Byte> nodeAseq = getVariantSequenceByNode(
					(FlexibleNode) selectedNode.getParent());
			ConcurrentHashMap<Integer, Byte> nodeBseq = getVariantSequenceByNode(selectedNode);

			Double[] selectedScores = new Double[3];
			selectedScores[0] = (double) get_num_leaves(selectedNode.getParent());
			selectedScores[1] = (double) get_num_leaves(selectedNode);
			ConcurrentHashMap<Integer, Byte> nodePseq = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq,
					selectedScores);

			String selectedBnid = (String) selectedNode.getAttribute(internalnode_nidname);

			if (selectedBnid == null) {
				selectedBnid = selectedNode.getTaxon().getId();
			}

			String selectAName = (String) selectedNode.getParent().getAttribute(internalnode_nidname);

			// Q-P pendent length local estimation
			double pqlen = 0.0;
			if (selectedScores[2] <= MinDoubleNumLimit)
				pqlen = 0.0;
			else /// selectedScores[2] > Double.MIN_VALUE
			{
				double scoreAB = computeNodeScore(nodeAseq, nodeBseq);
				FlexibleNode myNodeB = selectedNode;
				FlexibleNode myNodeA = myNodeB.getParent();
				/// iteratively consider upper branch of A¡¯s parent to A for scaling if
				/// selectedScores[2] > Double.MIN_VALUE and scoreAB <= MinDoubleNumLimit.
				while ((scoreAB <= MinDoubleNumLimit || myNodeB.getLength() <= MinDoubleNumLimit)
						&& !myNodeA.isRoot()) {
					myNodeB = myNodeA;
					myNodeA = myNodeB.getParent();
					scoreAB = computeNodeScore(getVariantSequenceByNode(myNodeA), getVariantSequenceByNode(myNodeB));
				}
				if (scoreAB > MinDoubleNumLimit && myNodeB.getLength() > MinDoubleNumLimit)
					pqlen = localEstimation(selectedScores[2], scoreAB, myNodeB.getLength());
				else {
					double p = selectedScores[2] / ((double) getAlignmentLength());
					pqlen = JC69(p);
				}
			}
            
			//modified on 20230523 by YYT to first consider A-P branch
			double original_branchAB = selectedNode.getLength();
			if (selectedScores[0] <= MinDoubleNumLimit) {
					ABQ_brlen[0] = 0.0;
					ABQ_brlen[1] = original_branchAB;
					ABQ_brlen[2] = pqlen;
			}
			else if (selectedScores[1] <= MinDoubleNumLimit) {
				ABQ_brlen[0] = original_branchAB;
				ABQ_brlen[1] = 0.0;
				ABQ_brlen[2] = pqlen;
			} else {
				double Pratio = selectedScores[0] / ((double) (selectedScores[0] + selectedScores[1]));
				double newNodePLength = original_branchAB * Pratio;
				double newNodeBLength = original_branchAB * (1.0 - Pratio);

				ABQ_brlen[0] = newNodePLength;
				ABQ_brlen[1] = newNodeBLength;
				ABQ_brlen[2] = pqlen;
			}

			if (printDisInfoOnScreen) {
				if (k == 0) // the best inserted branch
					System.out.println(qname + "\t" + "*" + selectAName + "-" + selectedBnid + "\t" + "ABQ_brlen: "
							+ ABQ_brlen[0] + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2]);
				else
					System.out.println(qname + "\t" + selectAName + "-" + selectedBnid + "\t" + "ABQ_brlen: "
							+ ABQ_brlen[0] + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2]);
			}

			/// for placement only
			if (otype.equals("placement")) {
				placementStrings += "[" + node2edge.get((FlexibleNode) mytree.getNode(selectedNodeBIndex)) + ", "
						+ ABQ_brlen[0] + ", " + ABQ_brlen[2] + "]";
				if (k < selectedNodeList.size() - 1)
					placementStrings += ",\n\t";
				else
					placementStrings += "\n\t";
			}

			if (k > 0)
				continue;

			// real add the query node to the Tree copy.
			MyFlexibleTree mynewTree = new MyFlexibleTree(mytree, true);
			this.copyAttributeFromOneNodeToAnother((FlexibleNode) mytree.getRoot(), (FlexibleNode) mynewTree.getRoot());

			mynewTree.beginTreeEdit();
			FlexibleNode selected_nodeB = (FlexibleNode) mynewTree.getNode(selectedNodeBIndex);

			FlexibleNode selected_nodeA = selected_nodeB.getParent();
			Taxon qtaxon = new Taxon(qname);
			FlexibleNode selected_nodeQ = new FlexibleNode(qtaxon);
			selected_nodeQ.setLength(pqlen);

			FlexibleNode selected_nodeP = new FlexibleNode();
			// set the attributes of newly added node.
			selected_nodeP.setAttribute(internalnode_nidname, pid);
			selected_nodeQ.setAttribute(internalnode_nidname, qname);
			if(selectedScores[0] <= MinDoubleNumLimit) // A-P is zero branch length, meaning that Q is
			{										   // inserted into A directly.
				selected_nodeP = selected_nodeA;
				selected_nodeA.addChild(selected_nodeQ);
			}
			else if (selectedScores[1] <= MinDoubleNumLimit) { // P-B is zero branch length, meaning that Q is inserted into
															// B directly.
				if (selected_nodeB.isExternal()) { // If B is leaf, cannot add the Q directly there, must add P node.
					selected_nodeP.setLength(original_branchAB);
					selected_nodeP.addChild(selected_nodeQ);
					selected_nodeA.addChild(selected_nodeP);

					selected_nodeB.setLength(0.0);
					selected_nodeA.removeChild(selected_nodeB);
					selected_nodeP.addChild(selected_nodeB);
				} else {
					selected_nodeB.addChild(selected_nodeQ);
				}
			} else {
				selected_nodeP.setLength(ABQ_brlen[0]);
				selected_nodeB.setLength(ABQ_brlen[1]);

				selected_nodeP.addChild(selected_nodeQ);
				selected_nodeA.addChild(selected_nodeP);

				selected_nodeA.removeChild(selected_nodeB);
				selected_nodeP.addChild(selected_nodeB);
			}

			mynewTree.endTreeEdit();

			// adding the query sequence for next insertion
			if (!seqIdxMap.containsKey(qname)) {
				seqIdxMap.put(qname, multationSequencesMap.size());
				multationSequencesMap.add(nodeQseq);
			}
			if (!seqIdxMap.containsKey(pid)) {
				seqIdxMap.put(pid, multationSequencesMap.size());
				multationSequencesMap.add(nodePseq);
			}

			// output the new added p node
			if (OUTPUT_PNODE) {
				writeMutation(pid, nodePseq, OUTPUT_FOLDER);
			}

			mynewTree.toAdoptNodes((FlexibleNode) mynewTree.getRoot());
			best_mynewTree = mynewTree;
		}
		/// to generate jplace file format info string
		if (otype.equals("placement")) {
			String placeInfo = "\t{\"p\":[\n\t" + placementStrings + "],\n\t" + "\"n\":[\"" + qname + "\"]}";
			placements[ii] = placeInfo;
			return mytree;
		}
		return best_mynewTree;
	}

	/// for fasta file
	public Tree addQuerySequence_global(String qname, String nodeQseq, String qid, String pid,
			boolean printDisInfoOnScreen, double[] ABQ_brlen, String otype, int ii) {
		// Travel all internal nodes for all possible nodeA-nodeB pairs

		minGlobalMemoryBranchScore = Double.MAX_VALUE;
		minGlobalMemoryBranchScoreNodeList.clear();
		/// if parallel computes branchscore is bigger than minGlobalMemoryBranchScore,
		/// the calculation will stop.
		ArrayList<Integer> nodeIdxAndScoreList = new ArrayList<Integer>(mytree.getNodeCount());
		//// initialize the nodeIdxAndScoreList to store nodeIdx
		for (int i = 0; i < mytree.getNodeCount(); ++i) {
			nodeIdxAndScoreList.add(i);
		}

		//// parallel to compute branchScore for every branch
		ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
			nodeIdxAndScoreList.parallelStream().forEach(nodeIdx -> {
				double score = computeBranchScore(nodeIdx, nodeQseq, qname);
			});
		});
		try {
			forkJoinTask.get();
		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		ArrayList<FlexibleNode> selectedNodeList = minGlobalMemoryBranchScoreNodeList;
		if (DEBUG)
			System.out.println(
					"minQScore: " + minGlobalMemoryBranchScore + ", selectedNodeList:" + selectedNodeList.size());

		//// try to remove ambiguous results
		selectedNodeList = reduceAmbiguousBranch(selectedNodeList, isMultiplePlacements);

		MyFlexibleTree best_mynewTree = null;

		String placementStrings = "";
		for (int k = 0; k < selectedNodeList.size(); ++k) {
			FlexibleNode selectedNode = selectedNodeList.get(k);

			int selectedNodeBIndex = selectedNode.getNumber();
			String nodeAseq = getStringSequenceByNode((FlexibleNode) selectedNode.getParent());
			String nodeBseq = getStringSequenceByNode(selectedNode);

			Double[] selectedScores = new Double[3];
			selectedScores[0] = (double) get_num_leaves(selectedNode.getParent());
			selectedScores[1] = (double) get_num_leaves(selectedNode);
			String nodePseq = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq, selectedScores);

			String selectedBnid = (String) selectedNode.getAttribute(internalnode_nidname);
			if (selectedBnid == null) {
				selectedBnid = selectedNode.getTaxon().getId();
			}
			String selectAName = (String) selectedNode.getParent().getAttribute(internalnode_nidname);

			// Q-P pendent length local estimation
			double pqlen = 0.0;
			if (selectedScores[2] <= MinDoubleNumLimit)
				pqlen = 0.0;
			else /// selectedScores[2] > Double.MIN_VALUE
			{
				double scoreAB = computeNodeScore(nodeAseq, nodeBseq);
				FlexibleNode myNodeB = selectedNode;
				FlexibleNode myNodeA = myNodeB.getParent();
				while ((scoreAB <= MinDoubleNumLimit || myNodeB.getLength() <= MinDoubleNumLimit)
						&& !myNodeA.isRoot()) {
					myNodeB = myNodeA;
					myNodeA = myNodeB.getParent();
					scoreAB = computeNodeScore(getStringSequenceByNode(myNodeA), getStringSequenceByNode(myNodeB));
				}
				if (scoreAB > MinDoubleNumLimit && myNodeB.getLength() > MinDoubleNumLimit)
					pqlen = localEstimation(selectedScores[2], scoreAB, myNodeB.getLength());
				else {
					double p = selectedScores[2] / ((double) getAlignmentLength());
					pqlen = JC69(p);
				}
			}

			double original_branchAB = selectedNode.getLength();
			if (selectedScores[0] <= MinDoubleNumLimit) {
				ABQ_brlen[0] = 0.0;
				ABQ_brlen[1] = original_branchAB;
				ABQ_brlen[2] = pqlen;
			}
			else if (selectedScores[1] <= MinDoubleNumLimit) {
			 
				ABQ_brlen[0] = original_branchAB;
				ABQ_brlen[1] = 0.0;
				ABQ_brlen[2] = pqlen;
			 
			} else {
				double Pratio = selectedScores[0] / ((double) (selectedScores[0] + selectedScores[1]));
				double newNodePLength = original_branchAB * Pratio;
				double newNodeBLength = original_branchAB * (1.0 - Pratio);

				ABQ_brlen[0] = newNodePLength;
				ABQ_brlen[1] = newNodeBLength;
				ABQ_brlen[2] = pqlen;
			}

			if (printDisInfoOnScreen) {
				if (k == 0)
					System.out.println(qname + "\t" + "*" + selectAName + "-" + selectedBnid + "\t" + "ABQ_brlen: "
							+ ABQ_brlen[0] + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2]);
				else
					System.out.println(qname + "\t" + selectAName + "-" + selectedBnid + "\t" + "ABQ_brlen: "
							+ ABQ_brlen[0] + "\t" + ABQ_brlen[1] + "\t" + ABQ_brlen[2]);
			}

			if (otype.equals("placement")) {
				placementStrings += "[" + node2edge.get((FlexibleNode) mytree.getNode(selectedNodeBIndex)) + ", "
						+ ABQ_brlen[0] + ", " + ABQ_brlen[2] + "]";
				if (k < selectedNodeList.size() - 1)
					placementStrings += ",\n\t";
				else
					placementStrings += "\n\t";
			}

			if (otype.equals("graft") || otype.equals("assembly")) {
				if (!node2placementSeqName.containsKey(selectedBnid)) {
					ArrayList<String> placementSeqName = new ArrayList<String>();
					node2placementSeqName.put(selectedBnid, placementSeqName);
				}
				node2placementSeqName.get(selectedBnid).add(qname);
			}

			if (k > 0)
				continue;

			// add the query node to the Tree copy.
			MyFlexibleTree mynewTree = new MyFlexibleTree(mytree, true);
			this.copyAttributeFromOneNodeToAnother((FlexibleNode) mytree.getRoot(), (FlexibleNode) mynewTree.getRoot());

			mynewTree.beginTreeEdit();
			FlexibleNode selected_nodeB = (FlexibleNode) mynewTree.getNode(selectedNodeBIndex);
			FlexibleNode selected_nodeA = selected_nodeB.getParent();
			Taxon qtaxon = new Taxon(qname);
			FlexibleNode selected_nodeQ = new FlexibleNode(qtaxon);
			selected_nodeQ.setLength(pqlen);

			FlexibleNode selected_nodeP = new FlexibleNode();
			// set the attributes of newly added node.
			selected_nodeP.setAttribute(this.internalnode_nidname, pid);
			selected_nodeQ.setAttribute(this.internalnode_nidname, qname);
			//double original_scoreAB = computeNodeScore(nodeAseq, nodeBseq);
			if (selectedScores[0] <= MinDoubleNumLimit) { // A-P is zero branch length, meaning that Q is
															// inserted into A directly.
				selected_nodeP = selected_nodeA;
				selected_nodeA.addChild(selected_nodeQ);
			}
			else if (selectedScores[1] <= MinDoubleNumLimit) { // P-B is zero branch length, meaning that Q is inserted into
															// B directly.
				if (selected_nodeB.isExternal()) { // If B is leaf, cannot add the Q directly there, must add P node.
					selected_nodeP.setLength(original_branchAB);
					selected_nodeP.addChild(selected_nodeQ);
					selected_nodeA.addChild(selected_nodeP);

					selected_nodeB.setLength(0.0);
					selected_nodeA.removeChild(selected_nodeB);
					selected_nodeP.addChild(selected_nodeB);
				} else {
					  
					selected_nodeB.addChild(selected_nodeQ);
				 
				}								  
			} else {
				selected_nodeP.setLength(ABQ_brlen[0]);
				selected_nodeB.setLength(ABQ_brlen[1]);

				selected_nodeP.addChild(selected_nodeQ);
				selected_nodeA.addChild(selected_nodeP);

				selected_nodeA.removeChild(selected_nodeB);
				selected_nodeP.addChild(selected_nodeB);
			}
			mynewTree.endTreeEdit();
			if (!seqIdxMap.containsKey(qname)) {
				seqIdxMap.put(qname, stringSequencesList.size());
				stringSequencesList.add(nodeQseq);
			}
			if (!seqIdxMap.containsKey(pid)) {
				seqIdxMap.put(pid, stringSequencesList.size());
				stringSequencesList.add(nodePseq);
			}
			if (OUTPUT_PNODE) {
				writeFASTA(pid, nodePseq, OUTPUT_FOLDER);
			}
			mynewTree.toAdoptNodes((FlexibleNode) mynewTree.getRoot());
			best_mynewTree = mynewTree;
		}
		if (otype.equals("placement")) {
			String placeInfo = "\t{\"p\":[\n\t" + placementStrings + "],\n\t" + "\"n\":[\"" + qname + "\"]}";
			placements[ii] = placeInfo;
			return mytree;
		}
		return best_mynewTree;
	}
	
	/// for fasta file
	public Tree addQuerySequence2Tree(FlexibleNode selectedNode, String qname, String nodeQseq, String qid, String pid, boolean printDisInfoOnScreen) 
	{		
		MyFlexibleTree best_mynewTree = null;
		int selectedNodeBIndex = selectedNode.getNumber();
		String nodeAseq = getStringSequenceByNode((FlexibleNode) selectedNode.getParent());
		String nodeBseq = getStringSequenceByNode(selectedNode);	

		Double[] selectedScores = new Double[3];
		selectedScores[0] = (double) get_num_leaves(selectedNode.getParent());
		selectedScores[1] = (double) get_num_leaves(selectedNode);
		String nodePseq = getStringAndScoreFromNodeABQSeq(nodeAseq, nodeBseq, nodeQseq, selectedScores);	

		String selectedBnid = (String) selectedNode.getAttribute(internalnode_nidname);
		if (selectedBnid == null) {
			selectedBnid = selectedNode.getTaxon().getId();
		}
		String selectAName = (String) selectedNode.getParent().getAttribute(internalnode_nidname);	

		// Q-P pendent length local estimation
		double pqlen = 0.0;
		if (selectedScores[2] <= MinDoubleNumLimit)
			pqlen = 0.0;
		else /// selectedScores[2] > Double.MIN_VALUE
		{
			double scoreAB = computeNodeScore(nodeAseq, nodeBseq);
			FlexibleNode myNodeB = selectedNode;
			FlexibleNode myNodeA = myNodeB.getParent();
			while ((scoreAB <= MinDoubleNumLimit || myNodeB.getLength() <= MinDoubleNumLimit)
					&& !myNodeA.isRoot()) {
				myNodeB = myNodeA;
				myNodeA = myNodeB.getParent();
				scoreAB = computeNodeScore(getStringSequenceByNode(myNodeA), getStringSequenceByNode(myNodeB));
			}
			if (scoreAB > MinDoubleNumLimit && myNodeB.getLength() > MinDoubleNumLimit)
				pqlen = localEstimation(selectedScores[2], scoreAB, myNodeB.getLength());
			else {
				double p = selectedScores[2] / ((double) getAlignmentLength());
				pqlen = JC69(p);
			}
		}				

		double original_branchAB = selectedNode.getLength();
		double[] ABQ_brlen = {0,0,0};
		if (selectedScores[0] <= MinDoubleNumLimit) {
			ABQ_brlen[0] = 0.0;
			ABQ_brlen[1] = original_branchAB;
			ABQ_brlen[2] = pqlen;
		}
		else if (selectedScores[1] <= MinDoubleNumLimit) {
			ABQ_brlen[0] = original_branchAB;
			ABQ_brlen[1] = 0.0;
			ABQ_brlen[2] = pqlen;
		} else {
			double Pratio = selectedScores[0] / ((double) (selectedScores[0] + selectedScores[1]));
			double newNodePLength = original_branchAB * Pratio;
			double newNodeBLength = original_branchAB * (1.0 - Pratio);

			ABQ_brlen[0] = newNodePLength;
			ABQ_brlen[1] = newNodeBLength;
			ABQ_brlen[2] = pqlen;
		}

		// add the query node to the Tree copy.
		MyFlexibleTree mynewTree = new MyFlexibleTree(mytree, true);
		copyAttributeFromOneNodeToAnother((FlexibleNode) mytree.getRoot(), (FlexibleNode) mynewTree.getRoot());

		mynewTree.beginTreeEdit();
		FlexibleNode selected_nodeB = (FlexibleNode) mynewTree.getNode(selectedNodeBIndex);
		FlexibleNode selected_nodeA = selected_nodeB.getParent();
		Taxon qtaxon = new Taxon(qname);
		FlexibleNode selected_nodeQ = new FlexibleNode(qtaxon);
		selected_nodeQ.setLength(pqlen);				

		FlexibleNode selected_nodeP = new FlexibleNode();	
		// set the attributes of newly added node.
		selected_nodeP.setAttribute(internalnode_nidname, pid);
		selected_nodeQ.setAttribute(internalnode_nidname, qname);
		//double original_scoreAB = computeNodeScore(nodeAseq, nodeBseq);
		if (selectedScores[0] <= MinDoubleNumLimit) { // A-P is zero branch length, meaning that Q is
														// inserted into A directly.
			selected_nodeP = selected_nodeA;
			selected_nodeA.addChild(selected_nodeQ);
		}
		else if (selectedScores[1] <= MinDoubleNumLimit) { // P-B is zero branch length, meaning that Q is inserted into
														// B directly.
			if (selected_nodeB.isExternal()) { // If B is leaf, cannot add the Q directly there, must add P node.
				selected_nodeP.setLength(original_branchAB);
				selected_nodeP.addChild(selected_nodeQ);
				selected_nodeA.addChild(selected_nodeP);

				selected_nodeB.setLength(0.0);
				selected_nodeA.removeChild(selected_nodeB);
				selected_nodeP.addChild(selected_nodeB);
			} else {
				selected_nodeB.addChild(selected_nodeQ);
			}
		} else {
			selected_nodeP.setLength(ABQ_brlen[0]);
			selected_nodeB.setLength(ABQ_brlen[1]);

			selected_nodeP.addChild(selected_nodeQ);
			selected_nodeA.addChild(selected_nodeP);

			selected_nodeA.removeChild(selected_nodeB);
			selected_nodeP.addChild(selected_nodeB);
		}
		mynewTree.endTreeEdit();
		
		if (!seqIdxMap.containsKey(qname)) {
			seqIdxMap.put(qname, stringSequencesList.size());
			stringSequencesList.add(nodeQseq);
		}
		if (!seqIdxMap.containsKey(pid)) {
			seqIdxMap.put(pid, stringSequencesList.size());
			stringSequencesList.add(nodePseq);
		}
		if (OUTPUT_PNODE) {
			writeFASTA(pid, nodePseq, OUTPUT_FOLDER);
		}
		mynewTree.toAdoptNodes((FlexibleNode) mynewTree.getRoot());
		best_mynewTree = mynewTree;	
	
		return best_mynewTree;
	}
	
	/// compute the branchscore that nodeQseq difference from both nodeA and nodeB

	//// for fasta file
	public double computeBranchScore(int nodeIdx, String nodeQseq, String qname) {
		FlexibleNode nodeB = (FlexibleNode) mytree.getNode(nodeIdx);
		return computeBranchScore(nodeB, nodeQseq, qname);
	}

	public double computeBranchScore(FlexibleNode nodeB, String nodeQseq, String qname) {
		// We use this annotation:
		// nodeA
		// |
		// nodeP
		// | \
		// | \
		// nodeB nodeQ

		if (nodeB.isRoot())
			return Double.MAX_VALUE;

		// nodeA-nodeB pairs
		FlexibleNode nodeA = (FlexibleNode) nodeB.getParent();
		if (nodeA == null) {
			System.out.println("nodeA is null, and its nodeB is " + node2seqName.get(nodeB));
		}

		String nodeAseq = getStringSequenceByNode(nodeA);
		String nodeBseq = getStringSequenceByNode(nodeB);

		if (nodeAseq == null) {
			if (DEBUG)
				System.out.println("nodeAseq cannot find seq - " + node2seqName.get(nodeA));

		} else if (nodeBseq == null) {
			if (DEBUG)
				System.out.println("nodeBseq cannot find seq - " + node2seqName.get(nodeB));
		}

		double score = 0;
		for (int i = 0; i < getAlignmentLength(); i++) {
			char a_i = nodeAseq.charAt(i);
			char b_i = nodeBseq.charAt(i);
			char c_i = nodeQseq.charAt(i);

			// get the substitution scores
			double score_ab = _used_scoreTable[(byte) a_i][(byte) b_i];
			double score_ac = _used_scoreTable[(byte) a_i][(byte) c_i];
			double score_bc = _used_scoreTable[(byte) b_i][(byte) c_i];

			if (a_i != b_i && a_i != c_i && b_i != c_i) { // ATC
				score += (score_ac + score_bc) / 2;
			} else if (a_i == b_i && b_i != c_i) { // AAT
				score += score_bc;
			}

			/// stop when score larger than current minGlobalMemoryBranchScore
			if (score - minGlobalMemoryBranchScore > MinDoubleNumLimit)
				return score;
		}

		lock.lock();
		if (minGlobalMemoryBranchScore - score > MinDoubleNumLimit) {
			minGlobalMemoryBranchScore = score;
			minGlobalMemoryBranchScoreNodeList.clear();
			minGlobalMemoryBranchScoreNodeList.add(nodeB);
		} else if (Math.abs(minGlobalMemoryBranchScore - score) <= MinDoubleNumLimit) {
			minGlobalMemoryBranchScoreNodeList.add(nodeB);
		}
		lock.unlock();

		return score;
	}

	public double computeBranchScore(String nodeAseq, String nodeBseq, String nodeQseq) {
		double score = 0;
		for (int i = 0; i < getAlignmentLength(); i++) {
			char a_i = nodeAseq.charAt(i);
			char b_i = nodeBseq.charAt(i);
			char c_i = nodeQseq.charAt(i);

			double score_ab = _used_scoreTable[(byte) a_i][(byte) b_i];
			double score_ac = _used_scoreTable[(byte) a_i][(byte) c_i];
			double score_bc = _used_scoreTable[(byte) b_i][(byte) c_i];

			if (a_i != b_i && a_i != c_i && b_i != c_i) { // ATC
				score += (score_ac + score_bc) / 2;
			} else if (a_i == b_i && b_i != c_i) { // AAT
				score += score_bc;
			}
		}

		return score;
	}

	//// compute the branch score and generate the P node string.
	public String getStringAndScoreFromNodeABQSeq(String aSeq, String bSeq, String cSeq, Double[] scores) {
		scores[0] = 0.0;
		scores[1] = 0.0;
		scores[2] = 0.0;
		StringBuilder pSeq = new StringBuilder(aSeq);

		// scores[3] is the difference between A and B
		for (int i = 0; i < getAlignmentLength(); i++) {
			char a_i = aSeq.charAt(i);
			char b_i = bSeq.charAt(i);
			char c_i = cSeq.charAt(i);

			double score_ac = _used_scoreTable[(byte) a_i][(byte) c_i];
			double score_bc = _used_scoreTable[(byte) b_i][(byte) c_i];
			double score_ab = _used_scoreTable[(byte) a_i][(byte) b_i];

			if (a_i == b_i && a_i == c_i) {
				continue;
				// do nothing
			} else if (a_i != b_i && a_i != c_i && b_i != c_i) { // ATC
				scores[2] += (score_ac + score_bc) / 2.0;
				if (score_ac > score_bc)
					pSeq.setCharAt(i, b_i);
			} else if (a_i == b_i && b_i != c_i) { // AAT
				scores[2] += score_bc;
				pSeq.setCharAt(i, a_i);
			} else if (a_i != b_i && b_i == c_i) { // ATT
				scores[0] += score_ab;
				pSeq.setCharAt(i, b_i);
			} else if (a_i == c_i && b_i != c_i) { // TAT
				scores[1] += score_ab;
				pSeq.setCharAt(i, a_i);
			} else if (DEBUG) {
				System.out.println("Unmatched ABQ type");
			}
		}
		return pSeq.toString();
	}

	//// for vcf file
	public double computeBranchScore(int nodeIdx, ConcurrentHashMap<Integer, Byte> nodeQseq, String qname) {
		FlexibleNode nodeB = (FlexibleNode) mytree.getNode(nodeIdx);
		return computeBranchScore(nodeB, nodeQseq, qname);
	}

	public double computeBranchScore(FlexibleNode nodeB, ConcurrentHashMap<Integer, Byte> nodeQseq, String qname) {
		// We use this annotation:
		// nodeA
		// |
		// nodeP
		// | \
		// | \
		// nodeB nodeQ

		if (nodeB.isRoot())
			return Double.MAX_VALUE;

		// nodeA-nodeB pairs
		FlexibleNode nodeA = (FlexibleNode) nodeB.getParent();
		ConcurrentHashMap<Integer, Byte> nodeAseq = getVariantSequenceByNode(nodeA);
		ConcurrentHashMap<Integer, Byte> nodeBseq = getVariantSequenceByNode(nodeB);

		double score = 0;
		HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
		mergeIdxSet.addAll(nodeAseq.keySet());
		mergeIdxSet.addAll(nodeBseq.keySet());

		for (Integer key : mergeIdxSet) {
			byte a = nodeAseq.containsKey(key) ? nodeAseq.get(key) : ref_sequence[key];
			byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
			byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

			if (a != b && a != c && b != c) { // ATC
				score += (_used_scoreTable[a][c] + _used_scoreTable[b][c]) / 2;
			} else if (a == b && b != c) { // AAT
				score += _used_scoreTable[b][c];
			}

			if (score - minGlobalMemoryBranchScore > MinDoubleNumLimit)
				return score;
		}

		lock.lock();
		if (minGlobalMemoryBranchScore - score > MinDoubleNumLimit) {
			minGlobalMemoryBranchScore = score;
			minGlobalMemoryBranchScoreNodeList.clear();
			minGlobalMemoryBranchScoreNodeList.add(nodeB);
		} else if (Math.abs(minGlobalMemoryBranchScore - score) <= MinDoubleNumLimit) {
			minGlobalMemoryBranchScoreNodeList.add(nodeB);
		}
		lock.unlock();

		return score;
	}

	public double computeBranchScore(ConcurrentHashMap<Integer, Byte> nodeAseq,
			ConcurrentHashMap<Integer, Byte> nodeBseq, ConcurrentHashMap<Integer, Byte> nodeQseq) {
		double score = 0;
		HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
		mergeIdxSet.addAll(nodeAseq.keySet());
		mergeIdxSet.addAll(nodeBseq.keySet());

		for (Integer key : mergeIdxSet) {
			byte a = nodeAseq.containsKey(key) ? nodeAseq.get(key) : ref_sequence[key];
			byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
			byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

			double score_ab = _used_scoreTable[a][b];
			double score_ac = _used_scoreTable[a][c];
			double score_bc = _used_scoreTable[b][c];

			if (a != b && a != c && b != c) { // ATC
				score += (score_ac + score_bc) / 2;
			} else if (a == b && b != c) { // AAT
				score += score_bc;
			}
		}
		return score;
	}

	public ConcurrentHashMap<Integer, Byte> getStringAndScoreFromNodeABQSeq(ConcurrentHashMap<Integer, Byte> nodeAseq,
			ConcurrentHashMap<Integer, Byte> nodeBseq, ConcurrentHashMap<Integer, Byte> nodeQseq, Double[] scores) {
		scores[0] = 0.0;
		scores[1] = 0.0;
		scores[2] = 0.0;
		ConcurrentHashMap<Integer, Byte> pSeq = new ConcurrentHashMap<Integer, Byte>();

		HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
		mergeIdxSet.addAll(nodeAseq.keySet());
		mergeIdxSet.addAll(nodeBseq.keySet());

		for (Integer key : mergeIdxSet) {
			byte a = nodeAseq.containsKey(key) ? nodeAseq.get(key) : ref_sequence[key];
			byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
			byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

			double score_ab = _used_scoreTable[(byte) a][(byte) b];
			double score_ac = _used_scoreTable[(byte) a][(byte) c];
			double score_bc = _used_scoreTable[(byte) b][(byte) c];

			byte placeCharacter = a; /// default place a
			if (a == b && a == c) {
				continue;
			} else if (a != b && a != c && b != c) { // ATC
				scores[2] += (score_ac + score_bc) / 2.0;
				if (score_ac > score_bc)
					placeCharacter = b;
			} else if (a == b && b != c) { // AAT
				scores[2] += score_bc;
			} else if (a != b && b == c) { // ATT
				scores[0] += score_ab;
				placeCharacter = b;
			} else if (a == c && b != c) { // TAT
				scores[1] += score_ab;
			}

			if (placeCharacter != ref_sequence[key]) {
				pSeq.put(key, placeCharacter);
			}
		}

		return pSeq;
	}

	public double computeNodeScore(FlexibleNode nodeB, String nodeQseq) {
		if (nodeB.isRoot())
			return Double.MAX_VALUE;

		String nodeBseq = getStringSequenceByNode(nodeB);
		double scoreNode = 0;
		for (int i = 0; i < getAlignmentLength(); i++) {
			char a_i = nodeBseq.charAt(i);
			char b_i = nodeQseq.charAt(i);
			scoreNode += _used_scoreTable[(byte) a_i][(byte) b_i];
		}

		return scoreNode;
	}

	public double computeNodeScore(String nodeBseq, String nodeQseq) {
		double score = 0;
		for (int i = 0; i < getAlignmentLength(); i++) {
			char a_i = nodeBseq.charAt(i);
			char b_i = nodeQseq.charAt(i);
			score += _used_scoreTable[(byte) a_i][(byte) b_i];
		}
		return score;
	}

	public double computeNodeScore(ConcurrentHashMap<Integer, Byte> nodeBseq,
			ConcurrentHashMap<Integer, Byte> nodeQseq) {
		double score = 0;
		HashSet<Integer> mergeIdxSet = new HashSet<Integer>(nodeQseq.keySet());
		mergeIdxSet.addAll(nodeBseq.keySet());

		for (Integer key : mergeIdxSet) {
			byte b = nodeBseq.containsKey(key) ? nodeBseq.get(key) : ref_sequence[key];
			byte c = nodeQseq.containsKey(key) ? nodeQseq.get(key) : ref_sequence[key];

			score += _used_scoreTable[b][c];
		}
		return score;
	}

	/// filter rules for multiple placements
	public ArrayList<FlexibleNode> reduceAmbiguousBranch(ArrayList<FlexibleNode> selectedNodeList) {
		return reduceAmbiguousBranch(selectedNodeList, false);
	}

	@SuppressWarnings("unchecked")
	public ArrayList<FlexibleNode> reduceAmbiguousBranch(ArrayList<FlexibleNode> selectedNodeList,
			Boolean returnMulti) {
		if (selectedNodeList.size() < 2)
			return selectedNodeList;
		ArrayList<FlexibleNode> finalSelectedNode = new ArrayList<FlexibleNode>();

		/// mapping slected nodes' parents to themselves
		HashMap<FlexibleNode, ArrayList<FlexibleNode>> parentCluster = new HashMap<FlexibleNode, ArrayList<FlexibleNode>>();
		for (int i = 0; i < selectedNodeList.size(); ++i) {
			FlexibleNode nodeB = selectedNodeList.get(i);
			FlexibleNode nodeA = (FlexibleNode) mytree.getParent(nodeB);
			if (!parentCluster.containsKey(nodeA)) {
				ArrayList<FlexibleNode> temp = new ArrayList<FlexibleNode>();
				temp.add(nodeB);
				parentCluster.put(nodeA, temp);
			} else {
				parentCluster.get(nodeA).add(nodeB);
			}
		}

		/// remove duplicate if potential placements are relationships of parent and
		/// child
		ArrayList<FlexibleNode> parentCluste_keySet = new ArrayList<FlexibleNode>(parentCluster.keySet());
		/// post transversal order
		parentCluste_keySet.sort((d1, d2) -> (d2.getNumber() - d1.getNumber()));
		for (int k = 0; k < parentCluste_keySet.size(); ++k) {
			FlexibleNode nodeA = parentCluste_keySet.get(k);
			if (parentCluster.containsKey(nodeA)) {
				ArrayList<FlexibleNode> temp = parentCluster.get(nodeA);
				for (int i = 0; i < temp.size(); ++i) {
					/// both parent node and its child node are potential placements
					/// select the node that contains more unique leaf nodes
					FlexibleNode nodeB = temp.get(i);
					if (parentCluster.containsKey(nodeB)) {
						int numLeavesA = get_num_leaves(nodeA);
						int numLeavesB = get_num_leaves(nodeB);
						if (numLeavesA < 2 * numLeavesB) {
							parentCluster.get(nodeA).remove(i);
							if (parentCluster.get(nodeA).size() < 1)
								parentCluster.remove(nodeA);
						} else
							parentCluster.remove(nodeB);
					}
				}
			}
		}

		/// select the node contain maximum children, if the same choose minimum height
		/// branch length
		/// for the case that several potential placements share the same parent
		for (FlexibleNode nodeA : parentCluster.keySet()) {
			ArrayList<FlexibleNode> temp = parentCluster.get(nodeA);
			int maxChildrenB = 0;
			ArrayList<Integer> selectIdxs = new ArrayList<Integer>();
			for (int i = 0; i < temp.size(); ++i) {
				if (temp.get(i).getChildCount() > maxChildrenB) {
					maxChildrenB = temp.get(i).getChildCount();
					selectIdxs.clear();
					selectIdxs.add(i);
				} else if (temp.get(i).getChildCount() == maxChildrenB) {
					selectIdxs.add(i);
				}
			}
			if (selectIdxs.size() < 2)
				finalSelectedNode.add(temp.get(selectIdxs.get(0)));
			else {
				double minHeight = Double.MAX_VALUE;
				for (int i = 0; i < selectIdxs.size(); ++i) {
					if (temp.get(selectIdxs.get(i)).getHeight() < minHeight) {
						minHeight = temp.get(selectIdxs.get(i)).getHeight();
					}
				}
				for (int i = 0; i < selectIdxs.size(); ++i) {
					if (Math.abs(minHeight - temp.get(selectIdxs.get(i)).getHeight()) <= MinDoubleNumLimit)
						finalSelectedNode.add(temp.get(selectIdxs.get(i)));
				}
			}
		}

		if (finalSelectedNode.size() < 2)
			return finalSelectedNode;
		selectedNodeList = (ArrayList<FlexibleNode>) finalSelectedNode.clone();
		finalSelectedNode.clear();

		/// return the best placement
		ArrayList<FlexibleNode> bestBranch = reduceBestBranch(selectedNodeList, parentCluster);
		if (returnMulti) {
			finalSelectedNode.add(bestBranch.get(0));
			// if multiple placements are allowed to print out, add all others
			for (int i = 0; i < selectedNodeList.size(); ++i) {
				if (selectedNodeList.get(i).getNumber() != bestBranch.get(0).getNumber()) {
					finalSelectedNode.add(selectedNodeList.get(i));
				}
			}
			return finalSelectedNode;
		} else
			return bestBranch;
	}

	/// filter different real multiple placements
	public ArrayList<FlexibleNode> reduceBestBranch(ArrayList<FlexibleNode> selectedNodeList,
			HashMap<FlexibleNode, ArrayList<FlexibleNode>> parentCluster) {
		if (selectedNodeList.size() < 2)
			return selectedNodeList;
		ArrayList<FlexibleNode> finalSelectedNode = new ArrayList<FlexibleNode>();

		//// select best answer ammong multiple placements
		/// select the maximum children then minimum height branch length
		int maxChildren = 0;
		for (int i = 0; i < selectedNodeList.size(); ++i) {
			FlexibleNode nodeB = selectedNodeList.get(i);
			FlexibleNode nodeA = (FlexibleNode) mytree.getParent(nodeB);

			int total = nodeA.getChildCount();
			if (parentCluster.containsKey(nodeB))
				total = total - 1; // unique children

			if (total > maxChildren) {
				maxChildren = total;
				finalSelectedNode.clear();
				finalSelectedNode.add(nodeB);
			} else if (maxChildren == total) {
				finalSelectedNode.add(nodeB);
			}
		}

		if (finalSelectedNode.size() < 2)
			return finalSelectedNode;
		selectedNodeList = new ArrayList<FlexibleNode>(finalSelectedNode);
		finalSelectedNode.clear();

		double minHeight = Double.MAX_VALUE;
		for (int i = 0; i < selectedNodeList.size(); ++i) {
			FlexibleNode nodeB = selectedNodeList.get(i);
			FlexibleNode nodeA = (FlexibleNode) mytree.getParent(nodeB);

			double total = nodeA.getHeight();

			if (minHeight - total > MinDoubleNumLimit) {
				minHeight = total;
				finalSelectedNode.clear();
				finalSelectedNode.add(nodeB);
			} else if (Math.abs(minHeight - total) <= MinDoubleNumLimit) {
				finalSelectedNode.add(nodeB);
			}
		}

		if (finalSelectedNode.size() > 1) {
			int randomIdx = (int) (Math.random() * finalSelectedNode.size());
			selectedNodeList = (ArrayList<FlexibleNode>) finalSelectedNode.clone();
			finalSelectedNode.clear();
			finalSelectedNode.add(selectedNodeList.get(randomIdx)); // random select
		}

		return finalSelectedNode;
	}

	public int get_num_leaves(FlexibleNode node) {
		if (!node.hasChildren())
			return 1;
		int num_leaves = 0;
		for (int i = 0; i < node.getChildCount(); ++i) {
			num_leaves += get_num_leaves(node.getChild(i));
		}
		return num_leaves;
	}

	public static void get_leaves(FlexibleNode node, HashSet<String> leaves) {
		if (!node.hasChildren()) {
			if (node2seqName.containsKey(node))
				leaves.add(node2seqName.get(node));
			else
				System.out.println("WARNNING: there is no sequence " + node.getNumber());
			return;
		}

		for (int i = 0; i < node.getChildCount(); ++i) {
			get_leaves(node.getChild(i), leaves);
		}
	}
	
	public static void get_leaves(FlexibleNode node, HashSet<String> leaves, HashMap<FlexibleNode, HashSet<String>> computedNodeLeaves) {
		if (!node.hasChildren()) {
			if (node2seqName.containsKey(node))
				leaves.add(node2seqName.get(node));
			else
				System.out.println("WARNNING: there is no sequence " + node.getNumber());
			return;
		}

		for (int i = 0; i < node.getChildCount(); ++i) {
			FlexibleNode mynode = node.getChild(i);
			if(computedNodeLeaves.containsKey(mynode)) leaves.addAll(computedNodeLeaves.get(mynode));
			else get_leaves(mynode, leaves);
		}
	}
	
	public static void get_leaves(FlexibleNode node, HashSet<String> leaves, ConcurrentHashMap<FlexibleNode, HashSet<String>> computedNodeLeaves) {
		if (!node.hasChildren()) {
			if (node2seqName.containsKey(node))
				leaves.add(node2seqName.get(node));
			else
				System.out.println("WARNNING: there is no sequence " + node.getNumber());
			return;
		}

		for (int i = 0; i < node.getChildCount(); ++i) {
			FlexibleNode mynode = node.getChild(i);
			if(computedNodeLeaves.containsKey(mynode)) leaves.addAll(computedNodeLeaves.get(mynode));
			else get_leaves(mynode, leaves);
		}
	}

	public static int get_leaves_exclude(FlexibleNode node, Set<FlexibleNode> excludeNodeList) {
		int count = 0;
		if (!excludeNodeList.contains(node))
		{
			if (node.hasChildren()) 
			{
				for (int i = 0; i < node.getChildCount(); ++i) 
				{
					count += get_leaves_exclude(node.getChild(i), excludeNodeList);
				}
			}
			else count = 1;
		}
		return count;
	}
	
	public static int get_leaves_exclude(FlexibleNode node, HashSet<FlexibleNode> samples_given, Set<FlexibleNode> excludeNodeList) {
		int count = 0;
		if (!excludeNodeList.contains(node))
		{
			if (node.hasChildren()) 
			{
				for (int i = 0; i < node.getChildCount(); ++i) 
				{
					count += get_leaves_exclude(node.getChild(i), excludeNodeList);
				}
			}
			else if(samples_given.contains(node)) count = 1;
		}
		return count;
	}
	

	public static int get_leaves_all(FlexibleNode node, HashMap<FlexibleNode, Integer> node2leaves) {
		int count = 0;
		if (node.hasChildren()) {
			for (int i = 0; i < node.getChildCount(); ++i) {
				count += get_leaves_all(node.getChild(i), node2leaves);
			}
		}
		else count = 1;
		node2leaves.put(node, count);
		return count;
	}
	
	public static int get_leaves_all(FlexibleNode node, HashSet<FlexibleNode> samples_given, HashMap<FlexibleNode, Integer> node2leaves) {
		int count = 0;
		if (node.hasChildren()) {
			for (int i = 0; i < node.getChildCount(); ++i) {
				count += get_leaves_all(node.getChild(i), node2leaves);
			}
		}
		else if(samples_given.contains(node)) count = 1;
	
		node2leaves.put(node, count);
		return count;
	}
	
	
	class MyFlexibleTree extends FlexibleTree {
		MyFlexibleTree(Tree t) {
			super(t);
		}

		MyFlexibleTree(Tree t, boolean keepAttribute) {
			super(t, keepAttribute);
		}

		/**
		 * Adopt a node hierarchy as its own. Only called by the
		 * FlexibleTree(FlexibleNode, TaxonList). This creates the node list and stores
		 * the nodes in post-traversal order.
		 */
		public void toAdoptNodes(FlexibleNode n) {
			super.adoptNodes(n);
		}

	}

	// JC69 model
	public static double JC69(double p) {
		double d = -3.0 * Math.log(1 - 4.0 * p / 3.0) / 4.0;
		return d;
	}

	// K2P model
	public static double K2P(String nodePseq, String nodeQseq) {
		int S = 0;
		int V = 0;
		int n = nodePseq.length();

		boolean precedingGapQ = true;

		for (int i = 0; i < n; i++) {
			if (nodeQseq.charAt(i) != '-') {
				precedingGapQ = false;
			}

			if (precedingGapQ && nodeQseq.charAt(i) == '-') {
				continue;
			}

			if (nodePseq.charAt(i) == nodeQseq.charAt(i)) {
				continue;
			} else if (nodePseq.charAt(i) == 'A' && nodeQseq.charAt(i) == 'G') {
				S++;
			} else if (nodePseq.charAt(i) == 'G' && nodeQseq.charAt(i) == 'A') {
				S++;
			} else if (nodePseq.charAt(i) == 'C' && nodeQseq.charAt(i) == 'T') {
				S++;
			} else if (nodePseq.charAt(i) == 'T' && nodeQseq.charAt(i) == 'C') {
				S++;
			} else {
				V++;
			}
		}

		double s = 1.0 * S / (1.0 * n);
		double v = 1.0 * V / (1.0 * n);
		double d = -1.0 * Math.log(1.0 - 2.0 * s - v) / 2.0 - 1.0 * Math.log(1.0 - 2.0 * v) / 4.0;
		// System.out.println(s + "\t" + v + "\t" + d);
		return d;
	}

	// local scale model
	public static double localEstimation(Integer[] selectedScores_int, double ABbranch) {
		if (selectedScores_int[0] + selectedScores_int[1] == 0) {
			return 0;
		}

		double d = ABbranch * ((double) selectedScores_int[2])
				/ ((double) (selectedScores_int[0] + selectedScores_int[1]));

		return d;
	}

	public static double localEstimation(Double[] selectedScores, double ABbranch) {
		if (selectedScores[0] + selectedScores[1] == 0) {
			return 0;
		}

		double d = ABbranch * selectedScores[2] / (selectedScores[0] + selectedScores[1]);

		return d;
	}

	public static double localEstimation(double score_c2ab, double score_ab, double ABbranch) {
		if (score_c2ab <= MinDoubleNumLimit) {
			return 0;
		}

		double d = ABbranch * score_c2ab / score_ab;

		return d;
	}

	// return common ancestor
	public static FlexibleNodeBranch returnCommonAncestor(FlexibleNode a, FlexibleNode b, Tree tree) {
		FlexibleNodeBranch nodeAndBranch = new TIPars().new FlexibleNodeBranch();
		if (a == b) {
			nodeAndBranch.a = a;
			nodeAndBranch.b = 0.0;
			return nodeAndBranch;
		} else {
			HashMap<FlexibleNode, Double> bAncestors = new HashMap<FlexibleNode, Double>();
			FlexibleNode parent = b;
			double dist = 0.0;
			bAncestors.put(b, dist);
			while (parent != null && !parent.isRoot()) {
				dist = dist + parent.getLength();
				bAncestors.put(parent.getParent(), dist);
				parent = parent.getParent();
			}
			parent = a;
			dist = 0.0;
			while (!bAncestors.containsKey(parent) && !tree.isRoot(parent)) {
				dist = dist + parent.getLength();
				parent = parent.getParent();
			}
			nodeAndBranch.a = parent;
			nodeAndBranch.b = dist;
			if (bAncestors.containsKey(parent))
				nodeAndBranch.b += bAncestors.get(parent);
			return nodeAndBranch;
		}
	}

	// remove single taxon
	public static FlexibleNodeBranch removeTaxon(FlexibleTree tree, FlexibleNode node) {
		try {
			MyFlexibleTree t = new TIPars().new MyFlexibleTree(tree, true);
			copyAttributeFromOneNodeToAnother((FlexibleNode) tree.getRoot(), (FlexibleNode) t.getRoot());

			t.beginTreeEdit();
			FlexibleNode n = (FlexibleNode) t.getNode(node.getNumber());
			FlexibleNode p = n.getParent();

			FlexibleNode pReturn = (FlexibleNode) tree.getNode(p.getNumber());
			FlexibleNode aReturn = (FlexibleNode) (pReturn.getParent());
			double true_dist1 = p.getLength() + n.getLength();
			double true_dist2 = n.getLength();

			p.removeChild(n);
			n.setParent(null);
			if (p.getChildCount() == 1) { // Remove this p node if if has less than 2 child nodes
				if (!p.isRoot()) {
					FlexibleNode b = (FlexibleNode) (p.getChild(0));
					FlexibleNode a = (FlexibleNode) (p.getParent());

					a.addChild(b);
					a.removeChild(p);
					b.setParent(a);
					double oldb2p = b.getLength();
					b.setLength(b.getLength() + p.getLength());
					p.setParent(null);
					t.endTreeEdit();
					FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(aReturn, true_dist1);
					result.tree = t;

					t.toAdoptNodes((FlexibleNode) t.getRoot());
					return result;
				} else {
					FlexibleNode b = (FlexibleNode) (p.getChild(0));
					for (int i = 0; i < b.getChildCount(); ++i) {
						FlexibleNode c = (FlexibleNode) (b.getChild(i));
						c.setLength(c.getLength() + b.getLength());
						p.addChild(c);
						c.setParent(p);
					}
					p.removeChild(b);
					b.setParent(null);
					t.endTreeEdit();
					FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(pReturn, true_dist1);
					result.tree = t;

					t.toAdoptNodes((FlexibleNode) t.getRoot());
					return result;
				}
			} else {
				t.endTreeEdit();
				FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(pReturn, true_dist2);
				result.tree = t;

				t.toAdoptNodes((FlexibleNode) t.getRoot());
				return result;// ? p, p.getLentth ?
			}
		} catch (Exception e) {
			e.printStackTrace();

			return null;
		}
	}

	// remove a list of taxa
	public static FlexibleNodeBranchList removeTaxon(FlexibleTree tree, ArrayList<FlexibleNode> nodes) {
		try {
			MyFlexibleTree t = new TIPars().new MyFlexibleTree(tree, true);
			copyAttributeFromOneNodeToAnother((FlexibleNode) tree.getRoot(), (FlexibleNode) t.getRoot());

			ArrayList<FlexibleNode> removeNodes = new ArrayList<FlexibleNode>();
			for (int i = 0; i < nodes.size(); ++i)
				removeNodes.add((FlexibleNode) t.getNode(nodes.get(i).getNumber()));

			FlexibleNodeBranchList flexibleNodeBranchList = new TIPars().new FlexibleNodeBranchList();
			flexibleNodeBranchList.nodeBranchList = new ArrayList<FlexibleNodeBranch>();

			t.beginTreeEdit();

			for (FlexibleNode node : removeNodes) {

				FlexibleNode n = node; // (FlexibleNode) t.getNode(node.getNumber());
				FlexibleNode p = n.getParent();
				System.out.println(p.getNumber());

				FlexibleNode pReturn = (FlexibleNode) tree.getNode(p.getNumber());
				String pName = (String) p.getAttribute(internalnode_nidname);
				String pReturnName = (String) pReturn.getAttribute(internalnode_nidname);
				if (!pName.equals(pReturnName)) {
					System.out.println("removeTaxon: " + pName + "|" + pReturnName);
					return null;
				}

				FlexibleNode aReturn = (FlexibleNode) (pReturn.getParent());
				double true_dist1 = p.getLength() + n.getLength();
				double true_dist2 = n.getLength();

				p.removeChild(n);
				n.setParent(null);
				if (p.getChildCount() == 1) { // Remove this p node if if has less than 2 child nodes
					FlexibleNode b = (FlexibleNode) (p.getChild(0));
					FlexibleNode a = (FlexibleNode) (p.getParent());

					a.addChild(b);
					b.setParent(a);
					double oldb2p = b.getLength();
					b.setLength(b.getLength() + p.getLength());

					a.removeChild(p);
					p.setParent(null);

					FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(aReturn, true_dist1);
					result.tree = null;
					flexibleNodeBranchList.nodeBranchList.add(result);
				} else {

					FlexibleNodeBranch result = new TIPars().new FlexibleNodeBranch(pReturn, true_dist2);
					result.tree = null;
					flexibleNodeBranchList.nodeBranchList.add(result);
				}
			}
			t.endTreeEdit();
			t.toAdoptNodes((FlexibleNode) t.getRoot());
			flexibleNodeBranchList.tree = t;
			return flexibleNodeBranchList;
		} catch (Exception e) {
			e.printStackTrace();

			return null;
		}
	}

	// remove a list of taxa and return a resulting tree
	public static MyFlexibleTree removeTaxonReturnTree(FlexibleTree tree, ArrayList<FlexibleNode> nodes) {
		try {
			MyFlexibleTree t = new TIPars().new MyFlexibleTree(tree, true);
			copyAttributeFromOneNodeToAnother((FlexibleNode) tree.getRoot(), (FlexibleNode) t.getRoot());

			ArrayList<FlexibleNode> removeNodes = new ArrayList<FlexibleNode>();
			for (int i = 0; i < nodes.size(); ++i) {
				removeNodes.add((FlexibleNode) t.getNode(nodes.get(i).getNumber()));
				String sequenceName1 = nodes.get(i).getTaxon().getId();
				String sequenceName2 = ((FlexibleNode) (t.getNode(nodes.get(i).getNumber()))).getTaxon().getId();
				if (!sequenceName1.equals(sequenceName2)) {
					System.out.println(sequenceName1 + " != " + sequenceName2);
					return null;
				}
				if (((FlexibleNode) (t.getNode(nodes.get(i).getNumber()))).isRoot()) {
					System.out.println(sequenceName2 + " original is root");
					return null;
				}
			}

			FlexibleNodeBranchList flexibleNodeBranchList = new TIPars().new FlexibleNodeBranchList();
			flexibleNodeBranchList.nodeBranchList = new ArrayList<FlexibleNodeBranch>();

			t.beginTreeEdit();
			for (FlexibleNode n : removeNodes) {
				if (n.isRoot()) {
					System.out.println(n.getTaxon().getId() + " is root");
					if (n.getParent() != null)
						System.out.println("root but not null");
					return null;
				}
				FlexibleNode p = n.getParent();
				if (p == null) {
					System.out.println(n.getTaxon().getId() + " no parent");
					return null;
				}
				if (p.getChildCount() > 2) {
					p.removeChild(n);
				} else if (p.getChildCount() == 2) {
					for (int j = 0; j < 2; ++j) {
						if (!p.getChild(j).equals(n)) {
							if (p.isRoot()) {
								//System.out.println(n.getTaxon().getId() + "'s parent is root and has only two children");
								//return null;
								//modified on 20211125
								FlexibleNode removeNode = p.getChild(j);
								for(int k=0; k<removeNode.getChildCount(); ++k)
								{
									FlexibleNode c = removeNode.getChild(k);
									double newbl = c.getLength() + removeNode.getLength();
									c.setLength(newbl);
									p.addChild(c);
									c.setParent(p);
								}
								p.removeChild(removeNode);
							}
							else
							{
								FlexibleNode a = (FlexibleNode) (p.getParent());
								FlexibleNode b = (FlexibleNode) (p.getChild(j));
								double newbl = b.getLength() + p.getLength();
								b.setLength(newbl);
								a.addChild(b);
								a.removeChild(p);
							}
						}
					}
				} else {
					System.out.println(n.getTaxon().getId() + " has no sister");
					return null;
				}
			}
			t.endTreeEdit();
			t.toAdoptNodes((FlexibleNode) t.getRoot());
			return t;
		} catch (Exception e) {
			e.printStackTrace();

			return null;
		}
	}

	// copy a node
	private static void copyAttributeFromOneNodeToAnother(FlexibleNode n, FlexibleNode m) {
		Iterator attnames = n.getAttributeNames();
		while (attnames != null && attnames.hasNext()) {
			String ahname = (String) attnames.next();
			m.setAttribute(ahname, n.getAttribute(ahname));
		}
	}

	// get string sequence
	private static String getStringSequenceByNode(FlexibleNode a) {
		if (node2seqName.containsKey(a)) {
			if (!seqIdxMap.containsKey((node2seqName.get(a)))) {
				System.out.println("seqIdxMap can not access " + (node2seqName.get(a)));
			}
			String seq = stringSequencesList.get(seqIdxMap.get((node2seqName.get(a))));
			return seq;
		} else {
			String selectedBnid = (String) a.getAttribute(internalnode_nidname);
			if (selectedBnid == null) {
				selectedBnid = a.getTaxon().getId();
			}
			System.out.println("node2seqName can not access " + a.getNumber() + "/" + selectedBnid);
			System.exit(-1);
			return null;
		}
	}

	private static String getStringSequenceByName(String a) {

		if (!seqIdxMap.containsKey(a)) {
			System.out.println("seqIdxMap can not access " + (a));
			System.exit(-1);
			return null;
		}
		String seq = stringSequencesList.get(seqIdxMap.get(a));
		return seq;
	}

	// get variants sequence
	private static ConcurrentHashMap<Integer, Byte> getVariantSequenceByNode(FlexibleNode a) {
		ConcurrentHashMap<Integer, Byte> seq = multationSequencesMap.get(seqIdxMap.get(node2seqName.get(a)));
		return seq;
	}
	
	private static ConcurrentHashMap<Integer, Byte> getVariantSequenceByName(String a) {
		ConcurrentHashMap<Integer, Byte> seq = multationSequencesMap.get(seqIdxMap.get(a));
		return seq;
	}

	public class FlexibleNodeBranch {
		public FlexibleNode a;
		public Double b;
		public Tree tree;
		public String nodeAName;

		public FlexibleNodeBranch(FlexibleNode a, Double b) {
			this.a = a;
			this.b = b;
		}

		public FlexibleNodeBranch() {
			this.a = null;
			this.b = null;
		}

		public void setNode(FlexibleNode a2) {
			this.a = a2;
		}

		public void setBrlen(Double b2) {
			this.b = b2;
		}

		public FlexibleNode getNode() {
			return this.a;
		}

		public Double getBrlen() {
			return this.b;
		}
	}

	public static int getAlignmentLength() {
		return sequence_character_length;
	}

	// read fasta file to a hashmap
	public static HashMap<Integer, String> readFastaFile2Alignment(String fn) {
		HashMap<Integer, String> alignmnetIdxList = new HashMap<Integer, String>();
		int startIndex = stringSequencesList.size();
		try {
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			String fasline = br2.readLine().trim();
			String desc;
			while (fasline != null) {
				StringBuffer sequence = new StringBuffer();
				if (fasline.matches("^>.+")) {
					desc = fasline;
					fasline = br2.readLine().trim();
					while (fasline != null && !fasline.matches("^>.+")) {
						sequence.append(fasline);
						fasline = br2.readLine();
					}
					desc = desc.replaceAll(">", "");
					String seqseq = sequence.toString();
					seqIdxMap.put(desc, startIndex);
					stringSequencesList.add(seqseq.toUpperCase());
					alignmnetIdxList.put(startIndex, desc);
					startIndex++;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.gc();
		return alignmnetIdxList;
	}

	// read fasta file
	public static void readFastaAlignmentFile(String fn) {
		int startIndex = stringSequencesList.size();
		try {
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			String fasline = br2.readLine().trim();
			String desc;
			while (fasline != null) {
				StringBuffer sequence = new StringBuffer();
				if (fasline.matches("^>.+")) {
					desc = fasline;
					fasline = br2.readLine().trim();
					while (fasline != null && !fasline.matches("^>.+")) {
						sequence.append(fasline);
						fasline = br2.readLine();
					}
					desc = desc.replaceAll(">", "");
					String seqseq = sequence.toString();
					if (sequence_character_length < 0)
						sequence_character_length = seqseq.length();
					seqIdxMap.put(desc, startIndex++);
					stringSequencesList.add(seqseq.toUpperCase());
				}
			}
		} catch (Exception e) {
			e.printStackTrace();

		}
		System.gc();
	}

	// read vcf file
	public static void readVCFAlignmentFile(String fn) {
		int startIndex = multationSequencesMap.size();
		boolean header_found = false;
		try {
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			//// read header
			while (!header_found) {
				String vcfline = br2.readLine().trim();
				if (vcfline == null || vcfline.isEmpty())
					continue;
				String[] words = vcfline.trim().split("\\s+");
				if (words.length > 2) {
					if (words[1].contains("POS")) {
						// Sample names start from the 10th word in the header
						for (int j = 9; j < words.length; j++) {
							seqIdxMap.put(words[j], startIndex + j - 9);
							multationSequencesMap.add(new ConcurrentHashMap<Integer, Byte>());
						}
						header_found = true;
					}
				}
			}

			if (!header_found) {
				System.out.println("Error! Incorrect VCF file. No Header");
				System.exit(-1);
			}
			int numberofThreads = THREAD_NUM; //Runtime.getRuntime().availableProcessors();
			// System.out.println("numberofThreads = " + numberofThreads);
			ArrayList<Integer> nodeIdxList = new ArrayList<Integer>(numberofThreads);
			for (int i = 0; i < numberofThreads; ++i)
				nodeIdxList.add(i);
			
			ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
				nodeIdxList.parallelStream().forEach(nodeIdx -> {
					String vcfline = null;
					try {
						while ((vcfline = br2.readLine()) != null) {
							if (vcfline.isEmpty())
								continue;
							String[] words = vcfline.trim().split("\\s+");
							if (words.length != 9 + multationSequencesMap.size() - startIndex) {
								System.out.println("Error! Incorrect VCF file.");
								System.exit(-1);
							}
							int variant_pos = Integer.parseInt(words[1]);
							byte ref_nuc = (byte) words[3].charAt(0);
							if (ref_sequence[variant_pos] == (byte) '\0') {
								ref_sequence[variant_pos] = ref_nuc;
							} else if (ref_nuc != ref_sequence[variant_pos]) {
								System.out.println("Error! Different referen sequence.");
								System.exit(-1);
							}

							String[] alleles = words[4].split(",");
							for (int j = 9; j < words.length; j++) {
								if (Character.isDigit(words[j].charAt(0))) {
									int allele_id = Integer.parseInt(words[j]);
									if (allele_id > 0) {
										byte allele = (byte) alleles[allele_id - 1].charAt(0);
										multationSequencesMap.get(startIndex + j - 9).put(variant_pos, allele);
									}
								}
							}
						}
					} catch (NumberFormatException | IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				});
			});
			forkJoinTask.get();
			br2.close();
			for (int i = 0; i < ref_sequence.length; ++i)
				if (ref_sequence[i] != (byte) '\0')
					sequence_character_length = i;
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.gc();
	}

	//// read vcf file to a hashmap
	public static HashMap<Integer, String> readVCFFile2Alignment(String fn) {
		HashMap<Integer, String> alignmnetIdxList = new HashMap<Integer, String>();
		int startIndex = multationSequencesMap.size();
		boolean header_found = false;
		try {
			BufferedReader br2 = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));

			//// read header
			while (!header_found) {
				String vcfline = br2.readLine().trim();
				if (vcfline == null || vcfline.isEmpty())
					continue;
				String[] words = vcfline.trim().split("\\s+");
				if (words.length > 2) {
					if (words[1].contains("POS")) {
						// Sample names start from the 10th word in the header
						for (int j = 9; j < words.length; j++) {
							seqIdxMap.put(words[j], startIndex + j - 9);
							multationSequencesMap.add(new ConcurrentHashMap<Integer, Byte>());
							alignmnetIdxList.put(startIndex + j - 9, words[j]);
						}
						header_found = true;
					}
				}
			}

			if (!header_found) {
				System.out.println("Error! Incorrect VCF file. No Header");
				System.exit(-1);
			}
			int numberofThreads = Runtime.getRuntime().availableProcessors();
			ArrayList<Integer> nodeIdxList = new ArrayList<Integer>(numberofThreads);
			for (int i = 0; i < numberofThreads; ++i)
				nodeIdxList.add(i);
			
			ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
				nodeIdxList.parallelStream().forEach(nodeIdx -> {
					String vcfline = null;
					try {
						while ((vcfline = br2.readLine()) != null) {
							if (vcfline == null)
								break;
							if (vcfline.isEmpty())
								continue;
							String[] words = vcfline.trim().split("\\s+");
							if (words.length > 1) {
								if (words.length != 9 + multationSequencesMap.size() - startIndex) {
									System.out.println("Error! Incorrect VCF file.");
									System.exit(-1);
								}
								int variant_pos = Integer.parseInt(words[1]);
								byte ref_nuc = (byte) words[3].charAt(0);
								if (ref_sequence[variant_pos] == (byte) '\0') {
									ref_sequence[variant_pos] = ref_nuc;
								} else if (ref_nuc != ref_sequence[variant_pos]) {
									System.out.println("Error! Different referen sequence.");
									System.exit(-1);
								}

								String[] alleles = words[4].split(",");
								for (int j = 9; j < words.length; j++) {
									if (Character.isDigit(words[j].charAt(0))) {
										int allele_id = Integer.parseInt(words[j]);
										if (allele_id > 0) {
											byte allele = (byte) alleles[allele_id - 1].charAt(0);
											multationSequencesMap.get(startIndex + j - 9).put(variant_pos, allele);
										}
									}
								}
							}
						}
					} catch (NumberFormatException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
				});
			});
			forkJoinTask.get();
			br2.close();
			for (int i = 0; i < ref_sequence.length; ++i)
				if (ref_sequence[i] != (byte) '\0')
					sequence_character_length = i;
		} catch (Exception e) {
			e.printStackTrace();
		}
		System.gc();
		return alignmnetIdxList;
	}

	public static String getFolder(String filename) {
		String dir = filename.substring(0, filename.lastIndexOf("/") + 1);
		return dir;
	}

	// write nwk tree or jplace file
	public static void writeToTree(Tree tree, String fn, String otype) {
		try {
			StringBuilder buffer = new StringBuilder();
			if (otype.equals("placement")) {
				buffer.append("{\n\t\"tree\": \"");
			}
			toNewick(tree, (FlexibleNode) tree.getRoot(), buffer, otype);
			buffer.append(';');

			if (otype.equals("placement")) {
				buffer.append("\",\n\t\"placements\": [\n");
				for (int i = 0; i < placements.length; i++) {
					buffer.append(placements[i]);
					if (i == placements.length - 1) {
						buffer.append("\n\t],\n");
					} else {
						buffer.append(",\n");
					}
				}
				buffer.append("\t\"metadata\": {\"info\": \"placement using TIPars\"},\n");
				buffer.append("\t\"version\": 201703,\n");
				buffer.append("\t\"fields\": [\"edge_num\", \"distal_length\", \"pendant_length\"\n\t]\n}\n");
			}

			PrintStream out = new PrintStream(new FileOutputStream(new File(fn)));
			out.println(buffer.toString());
			out.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void toNewick(Tree tree, FlexibleNode node, StringBuilder buffer, String otype) {
		if (tree.isExternal(node)) {
			String label = tree.getTaxonId(node.getNumber());
			buffer.append(label);
			appendLength(tree, node, buffer);
			appendEdgeNumber(buffer, node, otype);
		} else {
			buffer.append('(');
			int n = tree.getChildCount(node);
			for (int i = 0; i < n; i++) {
				toNewick(tree, (FlexibleNode) tree.getChild(node, i), buffer, otype);
				if (i == (n - 1)) {
					buffer.append(')');
					String label = (String) node.getAttribute(internalnode_nidname);
					if (node.isExternal())
						label = node.getTaxon().getId();
					if (label != null)
						buffer.append(label);
				} else {
					buffer.append(',');
				}
			}
			FlexibleNode parentNode = (FlexibleNode) tree.getParent(node);
			if (parentNode != null) {
				appendLength(tree, node, buffer);
				appendEdgeNumber(buffer, node, otype);
			}
		}
	}

	private static void appendEdgeNumber(StringBuilder buffer, FlexibleNode node, String otype) {
		if (otype.equals("placement")) {
			buffer.append('{');
			int edge = node2edge.get(node).intValue();
			buffer.append(edge);
			buffer.append('}');
		}
		;
	}

	private static void appendLength(Tree tree, FlexibleNode node, StringBuilder buffer) {
		if (tree.hasBranchLengths()) {
			buffer.append(':');
			buffer.append(tree.getBranchLength(node));
		}
	}

	// write nexus tree
	public static void writeNexusTree(Tree t, String fn) {
		try {
			PrintStream fw = new PrintStream(new FileOutputStream(new File(fn)));
			NexusExporter kne = new NexusExporter(fw); // export the tree *with attributes* to the output nexus file.
			Tree[] ts = new Tree[1];
			ts[0] = t;
			kne.exportTrees(ts, true);
			fw.close();
		} catch (Exception e) {
			e.printStackTrace();

		}
	}

	// write fasta
	public static void writeFASTA(String desc, String seq, String output_folder) {
		String output = output_folder + desc + ".fas";
		StringBuilder buffer = new StringBuilder();
		buffer.append(">");
		buffer.append(desc);
		buffer.append("\n");
		buffer.append(seq);
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
			out.println(buffer.toString());
			out.close();

		} catch (Exception e) {
			e.printStackTrace();

		}
	}

	public static void writeFASTAwithFullPath(String desc, String seq, String output_folder) {
		String output = output_folder;
		StringBuilder buffer = new StringBuilder();
		buffer.append(">");
		buffer.append(desc);
		buffer.append("\n");
		buffer.append(seq);
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
			out.println(buffer.toString());
			out.close();

		} catch (Exception e) {
			e.printStackTrace();

		}
	}

	public static void writeFASTA(ArrayList<String> desc, ArrayList<String> seq, String output_folder) {
		String output = output_folder;
		StringBuilder buffer = new StringBuilder();
		for (int i = 0; i < desc.size(); ++i) {
			buffer.append(">");
			buffer.append(desc.get(i));
			buffer.append("\n");
			buffer.append(seq.get(i));
			buffer.append("\n");
		}
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
			out.print(buffer.toString());
			out.close();

		} catch (Exception e) {
			e.printStackTrace();

		}
	}

	public static void writeFASTA(String[] desc, String[] seq, String output_folder) {
		String output = output_folder;
		StringBuilder buffer = new StringBuilder();
		for (int i = 0; i < desc.length; ++i) {
			buffer.append(">");
			buffer.append(desc[i]);
			buffer.append("\n");
			buffer.append(seq[i]);
			buffer.append("\n");
		}
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
			out.print(buffer.toString());
			out.close();

		} catch (Exception e) {
			e.printStackTrace();

		}
	}

	// write variants
	public static void writeMutation(String desc, ConcurrentHashMap<Integer, Byte> seqMutaion, String output_folder) {
		String output = output_folder + desc + ".vcf";
		StringBuilder buffer = new StringBuilder();
		buffer.append(desc);
		buffer.append("\n");
		for (ConcurrentHashMap.Entry<Integer, Byte> entry : seqMutaion.entrySet()) {
			System.out.println(entry.getKey());
			System.out.println(entry.getValue());
			buffer.append(ref_sequence[entry.getKey()] + "|" + entry.getKey() + "|" + entry.getValue() + "\n");
		}
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(output)));
			out.println(buffer.toString());
			out.close();

		} catch (Exception e) {
			e.printStackTrace();

		}
	}
	

	public class FlexibleNodeBranchList {
		public ArrayList<FlexibleNodeBranch> nodeBranchList = new ArrayList<FlexibleNodeBranch>();
		public Tree tree;
	}

	public static boolean isNumeric(String str) {
		for (int i = 0; i < str.length(); i++) {
			// System.out.println(str.charAt(i));
			if (!Character.isDigit(str.charAt(i))) {
				return false;
			}
		}
		return true;
	}

	// generate IUPAC nucleotide codes substitution table
	public static double[][] generate_nucleotide_nomenclature_scoreTable() {
		int size = 0;
		for (int i = 0; i < alphabet_nt.length; ++i)
			if ((byte) alphabet_nt[i] > size)
				size = (byte) alphabet_nt[i];
		size = size + 1;
		double[][] NodeScoreBigTable = new double[size][size];
		for (int i = 0; i < size; ++i)
			Arrays.fill(NodeScoreBigTable[0], 0);

		// {'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D', 'N',
		// '-'};

		for (int i = 0; i < alphabet_nt.length - 1; ++i) {
			for (int j = i + 1; j < alphabet_nt.length - 1; ++j) {
				int counts = _nucleotide_nomenclature.get((byte) alphabet_nt[i]).size()
						* _nucleotide_nomenclature.get((byte) alphabet_nt[j]).size();
				int commoncounts = 0;
				// if(counts == 1)
				// {
				for (byte key : _nucleotide_nomenclature.get((byte) alphabet_nt[i])) {
					if (_nucleotide_nomenclature.get((byte) alphabet_nt[j]).contains(key))
						commoncounts++;
				}
				if (commoncounts < 1)
					NodeScoreBigTable[alphabet_nt[i]][alphabet_nt[j]] = 1.0; // - commoncounts / (double)counts;
				// }
			}
		}

		for (int i = 0; i < alphabet_nt.length - 1; ++i) {
			for (int j = 0; j < i; ++j) {
				NodeScoreBigTable[alphabet_nt[i]][alphabet_nt[j]] = NodeScoreBigTable[alphabet_nt[j]][alphabet_nt[i]];
			}
		}
		
		//set gap penalty, no use
		//double gap_penalty = 1.0;
		//for (int i = 0; i < alphabet_nt.length - 1; ++i) 
		//	NodeScoreBigTable[alphabet_nt[alphabet_nt.length - 1]][alphabet_nt[i]] = gap_penalty;
		//for (int i = 0; i < alphabet_nt.length - 1; ++i) 
		//	NodeScoreBigTable[alphabet_nt[i]][alphabet_nt[alphabet_nt.length - 1]] = gap_penalty;

		return NodeScoreBigTable;
	}

	//// A: 0001
	//// T: 0010
	//// C: 0100
	//// G: 1000
	//// -: 0000
	public static BitSet encodeNucleotide(byte b_i, BitSet nucleotides) {
		if (!_nucleotide_nomenclature.containsKey(b_i))
			return nucleotides;
		for (byte key : _nucleotide_nomenclature.get((byte) b_i)) {
			switch (key) {
			case 'A':
				nucleotides.set(3);
				break;
			case 'T':
				nucleotides.set(2);
				break;
			case 'C':
				nucleotides.set(1);
				break;
			case 'G':
				nucleotides.set(0);
				break;
			default:
			}
		}
		return nucleotides;
	}

	public static BitSet encodeNucleotide(byte nucleotide) {
		BitSet nucleotides = new BitSet(4);
		if (!_nucleotide_nomenclature.containsKey(nucleotide))
			return nucleotides;
		for (byte key : _nucleotide_nomenclature.get(nucleotide)) {
			switch (key) {
			case 'A':
				nucleotides.set(3);
				break;
			case 'T':
				nucleotides.set(2);
				break;
			case 'C':
				nucleotides.set(1);
				break;
			case 'G':
				nucleotides.set(0);
				break;
			default:
			}
		}
		return nucleotides;
	}

	public static int NucleotideBitSet2Int(BitSet nucleotides) {
		int nucleotides_int = 0;
		for (int i = 0; i < 4; ++i) {
			if (nucleotides.get(i))
				nucleotides_int += Math.pow(2, i);
		}
		return nucleotides_int;
	}

	public static HashMap<BitSet, Byte> generate_nucleotide_nomenclature_characterTable() {
		HashMap<BitSet, Byte> nucleotide_nomenclature_map2char = new HashMap<BitSet, Byte>();
		Byte[] nucleotide_nomenclature_array = { 'A', 'C', 'G', 'T', 'R', 'Y', 'M', 'K', 'S', 'W', 'H', 'B', 'V', 'D',
				'N', '-' };
		for (int i = 0; i < nucleotide_nomenclature_array.length; ++i) {
			Byte nucleotide = nucleotide_nomenclature_array[i];
			BitSet bitSet = encodeNucleotide(nucleotide);
			nucleotide_nomenclature_map2char.put(bitSet, nucleotide);
		}
		return nucleotide_nomenclature_map2char;
	}

	// IUPAC nucleotide codes mapping
	public static HashMap<Byte, HashSet<Byte>> generate_nucleotide_nomenclature() {
		HashMap<Byte, HashSet<Byte>> nucleotide_nomenclature = new HashMap<Byte, HashSet<Byte>>();
		// A
		HashSet<Byte> tempA = new HashSet<Byte>();
		tempA.add((byte) 'A');
		nucleotide_nomenclature.put((byte) 'A', tempA);
		// C
		HashSet<Byte> tempC = new HashSet<Byte>();
		tempC.add((byte) 'C');
		nucleotide_nomenclature.put((byte) 'C', tempC);
		// G
		HashSet<Byte> tempG = new HashSet<Byte>();
		tempG.add((byte) 'G');
		nucleotide_nomenclature.put((byte) 'G', tempG);
		// T
		HashSet<Byte> tempT = new HashSet<Byte>();
		tempT.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'T', tempT);
		// R
		HashSet<Byte> tempR = new HashSet<Byte>();
		tempR.add((byte) 'A');
		tempR.add((byte) 'G');
		nucleotide_nomenclature.put((byte) 'R', tempR);
		// R
		HashSet<Byte> tempY = new HashSet<Byte>();
		tempY.add((byte) 'C');
		tempY.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'Y', tempY);
		// M
		HashSet<Byte> tempM = new HashSet<Byte>();
		tempM.add((byte) 'A');
		tempM.add((byte) 'C');
		nucleotide_nomenclature.put((byte) 'M', tempM);
		// K
		HashSet<Byte> tempK = new HashSet<Byte>();
		tempK.add((byte) 'G');
		tempK.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'K', tempK);
		// S
		HashSet<Byte> tempS = new HashSet<Byte>();
		tempS.add((byte) 'C');
		tempS.add((byte) 'G');
		nucleotide_nomenclature.put((byte) 'S', tempS);
		// W
		HashSet<Byte> tempW = new HashSet<Byte>();
		tempW.add((byte) 'A');
		tempW.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'W', tempW);
		// H
		HashSet<Byte> tempH = new HashSet<Byte>();
		tempH.add((byte) 'A');
		tempH.add((byte) 'C');
		tempH.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'H', tempH);
		// B
		HashSet<Byte> tempB = new HashSet<Byte>();
		tempB.add((byte) 'C');
		tempB.add((byte) 'G');
		tempB.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'B', tempB);
		// V
		HashSet<Byte> tempV = new HashSet<Byte>();
		tempV.add((byte) 'A');
		tempV.add((byte) 'C');
		tempV.add((byte) 'G');
		nucleotide_nomenclature.put((byte) 'V', tempV);
		// D
		HashSet<Byte> tempD = new HashSet<Byte>();
		tempD.add((byte) 'A');
		tempD.add((byte) 'G');
		tempD.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'D', tempD);
		// N
		HashSet<Byte> tempN = new HashSet<Byte>();
		tempN.add((byte) 'A');
		tempN.add((byte) 'C');
		tempN.add((byte) 'G');
		tempN.add((byte) 'T');
		nucleotide_nomenclature.put((byte) 'N', tempN);
		// -
		HashSet<Byte> tempGap = new HashSet<Byte>();
		tempGap.add((byte) '-');
		nucleotide_nomenclature.put((byte) '-', tempGap);

		return nucleotide_nomenclature;
	}	
	
	//generate blosum62 substitution table
    public static double[][] generate_aminoacid_scoreTable()
    {
    	int size = 0;
    	for(int i=0; i<alphabet_aa.length; ++i)
    		if((byte)alphabet_aa[i] > size) size = (byte)alphabet_aa[i]; //change char to integer based on ascii code
    	size = size + 1;
    	double[][] NodeScoreBigTable = new double[size][size];
    	for(int i=0; i<size; ++i) Arrays.fill(NodeScoreBigTable[0], 0);

    	double[][] blosum62 = {
    	    	//A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V  B  Z  X  -
    			{4,-1,-2,-2,0,-1,-1,0,-2,-1,-1,-1,-1,-2,-1,1,0,-3,-2,0,-2,-1,0,-4},
    			{-1,5,0,-2,-3,1,0,-2,0,-3,-2,2,-1,-3,-2,-1,-1,-3,-2,-3,-1,0,-1,-4},
    			{-2,0,6,1,-3,0,0,0,1,-3,-3,0,-2,-3,-2,1,0,-4,-2,-3,3,0,-1,-4},
    			{-2,-2,1,6,-3,0,2,-1,-1,-3,-4,-1,-3,-3,-1,0,-1,-4,-3,-3,4,1,-1,-4},
    			{0,-3,-3,-3,9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4},
    			{-1,1,0,0,-3,5,2,-2,0,-3,-2,1,0,-3,-1,0,-1,-2,-1,-2,0,3,-1,-4},
    			{-1,0,0,2,-4,2,5,-2,0,-3,-3,1,-2,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4},
    			{0,-2,0,-1,-3,-2,-2,6,-2,-4,-4,-2,-3,-3,-2,0,-2,-2,-3,-3,-1,-2,-1,-4},
    			{-2,0,1,-1,-3,0,0,-2,8,-3,-3,-1,-2,-1,-2,-1,-2,-2,2,-3,0,0,-1,-4},
    			{-1,-3,-3,-3,-1,-3,-3,-4,-3,4,2,-3,1,0,-3,-2,-1,-3,-1,3,-3,-3,-1,-4},
    			{-1,-2,-3,-4,-1,-2,-3,-4,-3,2,4,-2,2,0,-3,-2,-1,-2,-1,1,-4,-3,-1,-4},
    			{-1,2,0,-1,-3,1,1,-2,-1,-3,-2,5,-1,-3,-1,0,-1,-3,-2,-2,0,1,-1,-4},
    			{-1,-1,-2,-3,-1,0,-2,-3,-2,1,2,-1,5,0,-2,-1,-1,-1,-1,1,-3,-1,-1,-4},
    			{-2,-3,-3,-3,-2,-3,-3,-3,-1,0,0,-3,0,6,-4,-2,-2,1,3,-1,-3,-3,-1,-4},
    			{-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4,7,-1,-1,-4,-3,-2,-2,-1,-2,-4},
    			{1,-1,1,0,-1,0,0,0,-1,-2,-2,0,-1,-2,-1,4,1,-3,-2,-2,0,0,0,-4},
    			{0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,1,5,-2,-2,0,-1,-1,0,-4},
    			{-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1,1,-4,-3,-2,11,2,-3,-4,-3,-2,-4},
    			{-2,-2,-2,-3,-2,-1,-2,-3,2,-1,-1,-2,-1,3,-3,-2,-2,2,7,-1,-3,-2,-1,-4},
    			{0,-3,-3,-3,-1,-2,-2,-3,-3,3,1,-2,1,-1,-2,-2,0,-3,-1,4,-3,-2,-1,-4},
    			{-2,-1,3,4,-3,0,1,-1,0,-3,-4,0,-3,-3,-2,0,-1,-4,-3,-3,4,1,-1,-4},
    			{-1,0,0,1,-3,3,4,-2,0,-3,-3,1,-1,-3,-1,0,-1,-3,-2,-2,1,4,-1,-4},
    			{0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-2,0,0,-2,-1,-1,-1,-1,-1,-4},
    			{-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,1},   			
    	};
    	
    	for(int i=0; i<alphabet_aa.length; ++i)
    	{
    		for(int j=0; j<alphabet_aa.length; ++j)
    		{
	    		NodeScoreBigTable[alphabet_aa[i]][alphabet_aa[j]] = blosum62[i][j];
    		}
    	}
    	 
    	return NodeScoreBigTable;
    }
    
	private static void get_parents2root(FlexibleNode node, ArrayList<FlexibleNode> res) {
		res.add(node);
		if (!node.isRoot() && node.getParent() != null)
			get_parents2root(node.getParent(), res);
	}
	
	private static char getCharWithMaxPosibility(char ref_base, HashMap<Byte, Double> counts, double baseDepthRate) {
		double maxPosibility = 0;
		byte max_base = 0;
		//ArrayList<Byte> all_base = new ArrayList<Byte>();
		for (Byte base : counts.keySet()) {
			 double posibility = counts.get(base); 
			 if(posibility - maxPosibility > MinDoubleNumLimit) 
			 { 
				 max_base = base; 
				 maxPosibility = posibility;
			 }
			 //all_base.add(base);
		}
		
		if(maxPosibility < baseDepthRate)
		{
			return ref_base;
		}
		else
		{
			return (char) max_base;
		}
	}
	
	private static ArrayList<FlexibleNode> getCommonAncestor(ArrayList<String> readNames, HashMap<String, FlexibleNode> seqName2node)
	{
		HashMap<FlexibleNode, Integer> ancestorMap = new HashMap<FlexibleNode, Integer>();
		
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
		for(int i=0; i<readNames.size(); ++i)
		{
			FlexibleNode node = seqName2node.get(readNames.get(i));
			ancestorList.clear(); get_parents2root(node, ancestorList);
			for(int j=0; j<ancestorList.size(); ++j)
			{
				if(ancestorMap.containsKey(ancestorList.get(j)))
				{
					int count = ancestorMap.get(ancestorList.get(j));
					ancestorMap.replace(ancestorList.get(j), count + 1);
				}
				else
				{
					ancestorMap.put(ancestorList.get(j), 1);
				}
			}
		}

		return new ArrayList<FlexibleNode>(ancestorMap.keySet());
	}
	
	private static ConcurrentHashMap<FlexibleNode, Integer> getCommonAncestorAndCount(ArrayList<String> readNames, HashMap<String, FlexibleNode> seqName2node)
	{
		ConcurrentHashMap<FlexibleNode, Integer> ancestorMap = new ConcurrentHashMap<FlexibleNode, Integer>();
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
		for(int i=0; i<readNames.size(); ++i)
		{
			FlexibleNode currentnode = seqName2node.get(readNames.get(i));
			ancestorList.clear(); get_parents2root(currentnode, ancestorList);
			ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
				ancestorList.parallelStream().forEach(node -> 
				{
					if(ancestorMap.containsKey(node))
					{
						int count = ancestorMap.get(node);
						ancestorMap.replace(node, count + 1);
					}
					else
					{
						ancestorMap.put(node, 1);
					}
				});
			});
			try {
				forkJoinTask.get();
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}

		return ancestorMap;
	}
	

	
	private static HashMap<FlexibleNode, Integer> getCommonAncestorAndCount_exclude(ArrayList<String> readNames, HashMap<String, FlexibleNode> seqName2node, Set<FlexibleNode> excludeNodeList)
	{
		HashMap<FlexibleNode, Integer> ancestorMap = new HashMap<FlexibleNode, Integer>();
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
		for(int i=0; i<readNames.size(); ++i)
		{
			String seqname = readNames.get(i);
			FlexibleNode currentnode = seqName2node.get(seqname);
			ancestorList.clear(); get_parents2root(currentnode, ancestorList);
			//the order of ancestorList is from currentnode up to root
			//along the path, once if there is an assignment, the upper path could not add 1.
			for(int j=0; j<ancestorList.size(); ++j) //can not set to be parallel due to the order
			{
				FlexibleNode node = ancestorList.get(j);
				if(excludeNodeList.contains(node)) break;
				if(ancestorMap.containsKey(node))
				{
					int count = ancestorMap.get(node);
					ancestorMap.replace(node, count + 1);
				}
				else
				{
					ancestorMap.put(node, 1);
				}
			};
		};

		return ancestorMap;
	}
	
	private static HashMap<FlexibleNode, Integer> getCommonAncestorAndCount_exclude_parallel(ArrayList<String> readNames, HashMap<String, FlexibleNode> seqName2node, Set<FlexibleNode> excludeNodeList)
	{
		CopyOnWriteArrayList<ArrayList<FlexibleNode>> all_ancestorList = new CopyOnWriteArrayList<ArrayList<FlexibleNode>>();
		ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
			readNames.parallelStream().forEach(seqname  -> //to see how to make it parallel such that ancestorMap.replace(node, count + 1) is atomic operation
			{
				//String seqname = readNames.get(i);
				FlexibleNode currentnode = seqName2node.get(seqname);
				ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
				get_parents2root(currentnode, ancestorList);
				all_ancestorList.add(ancestorList);
			});
		});
		try {
			forkJoinTask.get();
		} catch (InterruptedException | ExecutionException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		HashMap<FlexibleNode, Integer> ancestorMap = new HashMap<FlexibleNode, Integer>();
		for(int i=0; i<all_ancestorList.size(); ++i)
		{
			ArrayList<FlexibleNode> ancestorList = all_ancestorList.get(i);
			//the order of ancestorList is from currentnode up to root
			//along the path, once if there is an assignment, the upper path could not add 1.
			for(int j=0; j<ancestorList.size(); ++j) //can not set to be parallel due to the order
			{
				FlexibleNode node = ancestorList.get(j);
				if(excludeNodeList.contains(node)) break;
				if(ancestorMap.containsKey(node))
				{
					ancestorMap.replace(node, ancestorMap.get(node)+1);
				}
				else
				{
					ancestorMap.put(node, 1);
				}
			};
		};

		return ancestorMap;
	}
	
	public static HashMap<String, ArrayList<String>> readCladeSample(String fn, HashSet<String> taxaNames) {
		HashMap<String, ArrayList<String>> clade2samples = new HashMap<String, ArrayList<String>>();

		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			String oneline = null;
			//// read header
			while ((oneline = br.readLine()) != null) {
				if (oneline.isEmpty())
					continue;
				String[] words = oneline.trim().split("\\s+");
				if (words.length == 2) {
					if (!clade2samples.containsKey(words[0]))
						clade2samples.put(words[0], new ArrayList<String>());

					if (taxaNames.contains(words[1]))
						clade2samples.get(words[0]).add(words[1]);
					else
						System.out.printf("WARNNING: sample of %s with lineage/clade $s is not in tree taxa", words[1], words[0]);
				} else {
					System.out.println("WARNNING: get a wrong format clade sample " + oneline);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return clade2samples;
	}

	public static void writeAssignClade(String fn, HashMap<String, String> assignmentMap) {
		StringBuilder buffer = new StringBuilder();
		for (HashMap.Entry<String, String> entry : assignmentMap.entrySet()) {
			buffer.append(entry.getKey() + "\t" + entry.getValue() + "\n");
		}
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(fn)));
			out.print(buffer.toString());
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void writeAssignCladeAndTaxa(String fn, HashMap<FlexibleNode, String> assignedNode2Clade, HashMap<FlexibleNode, ScoreAndStringList> assignedNode2Samples) {
		StringBuilder buffer = new StringBuilder();
		buffer.append("annotated_node\tannotated_node_precedor\tdist_to_root\tpangolineage\tF1score\tsamples\n");
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
		HashMap<FlexibleNode, FlexibleNode> assignedNode2parent = new HashMap<FlexibleNode, FlexibleNode>();
		for (HashMap.Entry<FlexibleNode, String> entry : assignedNode2Clade.entrySet()) {
			//get the preceding annotated node
			FlexibleNode annotated_node = entry.getKey();
			ancestorList.clear(); get_parents2root(annotated_node, ancestorList);
			int j = 0;
			FlexibleNode annotated_node_precedor = ancestorList.get(ancestorList.size()-1); //default to be root
			for(j=1; j<ancestorList.size(); ++j) //the first item in ancestorList is annotated_node itself
			{
				if(assignedNode2Clade.containsKey(ancestorList.get(j)))
				{
					annotated_node_precedor = ancestorList.get(j);
					break;
				}
			}
			assignedNode2parent.put(annotated_node, annotated_node_precedor);
			
			//the parent of root is set to be null
			buffer.append(node2seqName.get(annotated_node) + "\t" + (annotated_node_precedor != annotated_node? node2seqName.get(annotated_node_precedor) : "null"));
			//distance to root
			double dist_to_root = 0;
			for(j=0; j<ancestorList.size(); ++j)
			{
				dist_to_root += ancestorList.get(j).getLength();
			}
			buffer.append("\t" + dist_to_root); //distance to root
			buffer.append("\t" + entry.getValue()); //pangolineage
			buffer.append("\t" + assignedNode2Samples.get(annotated_node).F1Score); //F1Score
			ArrayList<String> sampleList = assignedNode2Samples.get(annotated_node).SampleList;
			buffer.append("\t" + sampleList.get(0));
			for(int i = 1; i < sampleList.size(); ++i)
				buffer.append("," + sampleList.get(i));
			buffer.append("\n");
		}
		//print out the root if the root has not been annotated
		FlexibleNode myroot = (FlexibleNode) mytree.getRoot();
		if(!assignedNode2Clade.containsKey(myroot))
		{
			buffer.append(node2seqName.get(myroot) + "\tNA\t0\tnull\tnull\tnull\n");
		}
		
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(fn)));
			out.print(buffer.toString());
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void writeSeqName(String fn, ArrayList<String> seqnames) {
		StringBuilder buffer = new StringBuilder();
		for (int i = 0; i < seqnames.size(); ++i) {
			buffer.append(seqnames.get(i) + "\n");
		}
		try {
			PrintStream out = new PrintStream(new FileOutputStream(new File(fn)));
			out.print(buffer.toString());
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static ArrayList<String> readSeqName(String fn) {
		ArrayList<String> seqnames = new ArrayList<String>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			String oneline = null;
			//// read header
			while ((oneline = br.readLine()) != null) {
				if (oneline.isEmpty())
					continue;
				seqnames.add(oneline.trim());
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return seqnames;
	}

	// function to sort hashmap by values
    public static HashMap<String, Integer> sortByValue(HashMap<String, Integer> hm)
    {
        // Create a list from elements of HashMap
        List<Map.Entry<String, Integer> > list =
               new LinkedList<Map.Entry<String, Integer> >(hm.entrySet());
 
        // Sort the list
        Collections.sort(list, new Comparator<Map.Entry<String, Integer> >() {
            public int compare(Map.Entry<String, Integer> o1,
                               Map.Entry<String, Integer> o2)
            {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });
         
        // put data from sorted list to hashmap
        HashMap<String, Integer> temp = new LinkedHashMap<String, Integer>();
        for (Map.Entry<String, Integer> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }
    
    public static HashMap<String, Double> sortByDoubleValue(ConcurrentHashMap<String, Double> hm)
    {
        // Create a list from elements of HashMap
        List<Map.Entry<String, Double> > list =
               new LinkedList<Map.Entry<String, Double> >(hm.entrySet());
 
        // Sort the list
        Collections.sort(list, new Comparator<Map.Entry<String, Double> >() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2)
            {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });
         
        // put data from sorted list to hashmap
        HashMap<String, Double> temp = new LinkedHashMap<String, Double>();
        for (Map.Entry<String, Double> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }
    
    public static HashMap<String, Double> sortByDoubleValue(HashMap<String, Double> hm)
    {
        // Create a list from elements of HashMap
        List<Map.Entry<String, Double> > list =
               new LinkedList<Map.Entry<String, Double> >(hm.entrySet());
 
        // Sort the list
        Collections.sort(list, new Comparator<Map.Entry<String, Double> >() {
            public int compare(Map.Entry<String, Double> o1,
                               Map.Entry<String, Double> o2)
            {
                return (o2.getValue()).compareTo(o1.getValue());
            }
        });
         
        // put data from sorted list to hashmap
        HashMap<String, Double> temp = new LinkedHashMap<String, Double>();
        for (Map.Entry<String, Double> aa : list) {
            temp.put(aa.getKey(), aa.getValue());
        }
        return temp;
    }

    public HashMap<String, Double> computeF1Scores(HashMap<String, ArrayList<String>> clade2samples, HashMap<String, FlexibleNode> seqName2node, HashMap<FlexibleNode, Integer> node2leaves, boolean printDisInfoOnScreen)
    {
    	HashMap<String, Double> clade2F1Score = new HashMap<String, Double>();
    	for(Map.Entry<String, ArrayList<String>> entry : clade2samples.entrySet()) {
            String cladeName = entry.getKey();
            ArrayList<String> samples = entry.getValue();
            ConcurrentHashMap<FlexibleNode, Integer> ancestor2Count = getCommonAncestorAndCount(samples, seqName2node);
            AtomicReference<Double> maxF1Score = new AtomicReference<Double>(0.0);
            ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
            	ancestor2Count.entrySet().parallelStream().forEach(entryNode -> {
                	FlexibleNode currentNode = entryNode.getKey();
                	int count = entryNode.getValue();
                	int leaves = node2leaves.get(currentNode);
        			double Precision = count / (double) leaves;
        			double Recall = count / (double) samples.size();
        			double F1Score = 2*Precision*Recall / (Precision + Recall);
        			
        			if(F1Score > maxF1Score.get())
        			{
        				maxF1Score.set(F1Score);
        			}
        			
    			});
            });
            try {
				forkJoinTask.get();
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			clade2F1Score.put(cladeName, maxF1Score.get());
			
			if(printDisInfoOnScreen)
			{
				System.out.println("processing : " + cladeName);	
			}
        }
		
    	return clade2F1Score;
    }
    
        
    //tree clade/lineage annotation
	public HashMap<String, String> treeAnnotation(HashMap<String, ArrayList<String>> clade2samples, double fscore_min, boolean printDisInfoOnScreen) {
		
//		//sort clade2samples based on each clade samples number, no use
//		HashMap<String, Integer> clade2sampleNum = new HashMap<String, Integer>();
//		for(HashMap.Entry<String, ArrayList<String>> entry : clade2samples.entrySet())
//		{
//			clade2sampleNum.put(entry.getKey(), entry.getValue().size());
//		}
//		clade2sampleNum = sortByValue(clade2sampleNum);
		
		//do assignment based on common ancestors
		//Filter1: only process > overlapThreshold common ancestors, no use
		//Filter2: select the ancestor with the most concentration
		HashMap<String, FlexibleNode> seqName2node = setupHashtableOfseqName2node(mytree);
		HashMap<FlexibleNode, Integer> node2leaves = new HashMap<FlexibleNode, Integer>();
		
		HashSet<FlexibleNode> samples_given = getSampleWithAnnotation(clade2samples, seqName2node);
    	int num_taxa = get_leaves_all((FlexibleNode) mytree.getRoot(), samples_given, node2leaves);
    	
    	//sort the best F1Scores for each lineage
		HashMap<String, Double> clade2F1Score = sortByDoubleValue(computeF1Scores(clade2samples, seqName2node, node2leaves, printDisInfoOnScreen));
		HashMap<FlexibleNode, String> assignmentNode2Clade = new HashMap<FlexibleNode, String>();
		HashMap<String, String> assignmentMap = new HashMap<String, String>();

		for (Map.Entry<String, Double> entry : clade2F1Score.entrySet()) {
			if(printDisInfoOnScreen) System.out.println("Key = " + entry.getKey() + ", Value = " + entry.getValue());

			String cladeName = entry.getKey();
            ArrayList<String> samples = clade2samples.get(cladeName);
            
            HashMap<FlexibleNode, Integer> ancestor2Count = getCommonAncestorAndCount_exclude(samples, seqName2node, assignmentNode2Clade.keySet());
            
            //HashMap<FlexibleNode,Integer> ancestor2Count = getCommonAncestorAndCount_exclude_parallel(samples, seqName2node, assignmentNode2Clade.keySet());
            double[] maxF1Score_selectPrecisionRecall = new double[] {0,0,0};
			FlexibleNode[] selectNode = new FlexibleNode[] {null};
			int[] trueCount = new int[] {0};
			
			ForkJoinTask<?> forkJoinTask = forkJoinPool.submit(() -> {
				ancestor2Count.entrySet().parallelStream().forEach(entryNode -> {
					FlexibleNode currentNode = entryNode.getKey();
					
					int leaves = get_leaves_exclude(currentNode, samples_given, assignmentNode2Clade.keySet()); 
					int count = entryNode.getValue();//.get();
					double Precision = count / (double) leaves;
					double Recall = count / (double) samples.size();
					double F1Score = 2 * Precision * Recall / (Precision + Recall);
	    			
					lock.lock();
					if(F1Score > maxF1Score_selectPrecisionRecall[0] && !assignmentNode2Clade.containsKey(currentNode))
					{
						maxF1Score_selectPrecisionRecall[0] = F1Score;
						maxF1Score_selectPrecisionRecall[1] = Precision;
						maxF1Score_selectPrecisionRecall[2] = Recall;
						selectNode[0] = currentNode;
						trueCount[0] = count;
					}
					lock.unlock();	
	            });
			});
			try {
				forkJoinTask.get();
			} catch (InterruptedException | ExecutionException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
            if (selectNode[0] != null && maxF1Score_selectPrecisionRecall[0] >= fscore_min) {
            	assignmentNode2Clade.put(selectNode[0], cladeName);
				String selectNodeName = node2seqName.get(selectNode[0]);
				assignmentMap.put(selectNodeName, cladeName);
				if(printDisInfoOnScreen)
					System.out.printf("MSG: Assigned %s to node %s as %.2f%% descendants are %.2f%% total given lineage/clade with F1Score %.2f%n", cladeName,
						selectNodeName, maxF1Score_selectPrecisionRecall[1] * 100, maxF1Score_selectPrecisionRecall[2] * 100, maxF1Score_selectPrecisionRecall[0] * 100);
			} else {
				System.out.printf("WARNNING: cannot assign %s because its common ancestors have been already assigned%n", cladeName);
			}
        }
		
		//annotationStatistics(clade2samples, assignmentNode2Clade, seqName2node, null, printDisInfoOnScreen);
		//System.gc();
		return assignmentMap;
	}
	
	public class ScoreAndStringList {
		public double F1Score;
		public ArrayList<String> SampleList = new ArrayList<String>();
	}
	
	public HashMap<FlexibleNode, ScoreAndStringList> annotationStatistics(HashMap<String, ArrayList<String>> clade2samples, HashMap<FlexibleNode, String> assignmentNode2Clade, HashMap<String, FlexibleNode> seqName2node, HashSet<String> taxaNames, boolean printDisInfoOnScreen, boolean output_CollapsedTree_ConsistentTree) {
		///statistic
		int truePossitive = 0;
		int assignedSamples = 0;
		int totalSamples = 0;
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
		
		HashMap<FlexibleNode, ScoreAndStringList> assignedNode2Samples = new HashMap<FlexibleNode, ScoreAndStringList>();
		HashMap<FlexibleNode, Integer[]> assignedNode2TrueFalsePositive = new HashMap<FlexibleNode, Integer[]>();
		
		for(HashMap.Entry<String, ArrayList<String>> entry : clade2samples.entrySet())
		{
			String clade = entry.getKey();
			ArrayList<String> samples = entry.getValue();
			totalSamples += samples.size();
	
			for(int i=0; i<samples.size(); ++i)
			{
				if(seqName2node.containsKey(samples.get(i)))
				{
					FlexibleNode sampleNode = seqName2node.get(samples.get(i));
					ancestorList.clear(); get_parents2root(sampleNode, ancestorList);
					int j = 0;
					for(j=0; j<ancestorList.size(); ++j)
					{
						if(assignmentNode2Clade.containsKey(ancestorList.get(j)))
						{
							String predictClade = assignmentNode2Clade.get(ancestorList.get(j));
							assignedSamples++;
							
							if(!assignedNode2Samples.containsKey(ancestorList.get(j)))
							{
								ScoreAndStringList scoreAndStringList = new TIPars().new ScoreAndStringList();
								assignedNode2Samples.put(ancestorList.get(j), scoreAndStringList);
								Integer[] truePositiveNegativeCount = {0, 0};
								assignedNode2TrueFalsePositive.put(ancestorList.get(j), truePositiveNegativeCount);
							}
							assignedNode2Samples.get(ancestorList.get(j)).SampleList.add(samples.get(i));
							
							if(!predictClade.equals(clade))
							{
								if(printDisInfoOnScreen)
									System.out.printf("sample %s in clade/lineage %s has been assigned to different %s clade/lineage%n", samples.get(i), clade, predictClade);
								assignedNode2TrueFalsePositive.get(ancestorList.get(j))[1]++;
							}
							else
							{	
								truePossitive++;
								assignedNode2TrueFalsePositive.get(ancestorList.get(j))[0]++;
							}
							break;
						}
					}
					if(j == ancestorList.size())
					{
						if(printDisInfoOnScreen)
							System.out.printf("sample %s in clade/lineage %s has not been assigned to any clade/lineage%n", samples.get(i), clade);
					}
				}
				else
				{
					if(printDisInfoOnScreen)
						System.out.printf("WARNNING: given tree does not contain sample %s in clade/lineage %s%n", samples.get(i), clade);
				}
			}
		}
		
		for(HashMap.Entry<FlexibleNode, String> entry : assignmentNode2Clade.entrySet())
		{
			FlexibleNode annotated_node = entry.getKey();
			String clade = entry.getValue();
			double Precision = assignedNode2TrueFalsePositive.get(annotated_node)[0] / (double) (assignedNode2TrueFalsePositive.get(annotated_node)[0] +assignedNode2TrueFalsePositive.get(annotated_node)[1]);
			double Recall = assignedNode2TrueFalsePositive.get(annotated_node)[0] / (double) clade2samples.get(clade).size();
			assignedNode2Samples.get(annotated_node).F1Score = 2 * Precision * Recall / (Precision + Recall);
		}
		
		System.out.printf("INFO: %d out of %d clade/lineages have been annotated.%n", assignmentNode2Clade.size(), clade2samples.size());
		System.out.printf("INFO: %d out of %d samples have been assigned with true predictions %d.%n", assignedSamples, totalSamples, truePossitive);	
		
		//get back all taxa samples in the tree but not in the given clade2samples by users
		//but these taxa are not included in the calculation of F1Score
		if(taxaNames != null)
		{
			for(HashMap.Entry<String, ArrayList<String>> entry : clade2samples.entrySet())
			{
				ArrayList<String> samples = entry.getValue();
				for(int i=0; i<samples.size(); ++i)
				{
					taxaNames.remove(samples.get(i));
				}
			}
			for (String taxonname : taxaNames)
			{
				if(seqName2node.containsKey(taxonname))
				{
					FlexibleNode sampleNode = seqName2node.get(taxonname);
					ancestorList.clear(); get_parents2root(sampleNode, ancestorList);
					int j = 0;
					for(j=0; j<ancestorList.size(); ++j)
					{
						if(assignmentNode2Clade.containsKey(ancestorList.get(j)))
						{
							String predictClade = assignmentNode2Clade.get(ancestorList.get(j));
							if(printDisInfoOnScreen)
								System.out.printf("sample %s without labeled clade/lineage has been assigned to %s clade/lineage%n", taxonname, predictClade);
							if(!assignedNode2Samples.containsKey(ancestorList.get(j)))
							{
								ScoreAndStringList scoreAndStringList = new TIPars().new ScoreAndStringList();
								assignedNode2Samples.put(ancestorList.get(j), scoreAndStringList);
							}
							assignedNode2Samples.get(ancestorList.get(j)).SampleList.add(taxonname);
							break;
						}
					}
					if(j == ancestorList.size())
					{
						if(printDisInfoOnScreen)
							System.out.printf("sample %s without labeled clade/lineage has not been assigned to any clade/lineage%n", taxonname);
					}
				}
				else
				{
					if(printDisInfoOnScreen)
						System.out.printf("WARNNING: given tree does not contain sample %s without labeled clade/lineage%n", taxonname);
				}
			}
		}
		
		//write out the consistent label tree and inconsistent samples
		if(output_CollapsedTree_ConsistentTree)
			getConsistentLabelTree(clade2samples, printDisInfoOnScreen, assignmentNode2Clade);
		
		//write out the annotation collapsed tree
		if(output_CollapsedTree_ConsistentTree)
			getCollapsedTree(assignmentNode2Clade);
		
		return assignedNode2Samples;	
	}
	
	//get all sample with given annotation
	public HashSet<FlexibleNode> getSampleWithAnnotation(HashMap<String, ArrayList<String>> clade2samples, HashMap<String, FlexibleNode> seqName2node)
	{
		HashSet<FlexibleNode> samples_given = new HashSet<FlexibleNode>();
		for(HashMap.Entry<String, ArrayList<String>> entry : clade2samples.entrySet())
		{
			ArrayList<String> samples = entry.getValue();
			for(int i=0; i<samples.size(); ++i)
			{
				samples_given.add(seqName2node.get(samples.get(i)));
			}
		}
		return samples_given;
	}
	
	//remove any inconsistent or unassigned lineage/clade taxa from the tree 
	public Tree removeInconsistentSample(HashMap<String, ArrayList<String>> clade2samples, HashMap<String, String> assignmentMap, Tree tree, HashMap<String, String> removedSampleAndClades)
	{
		HashMap<String, FlexibleNode> seqName2node = setupHashtableOfseqName2node(tree);
		HashMap<FlexibleNode, String> assignmentNode2Clade = new HashMap<FlexibleNode, String>();
		
		for(HashMap.Entry<String,String> entry : assignmentMap.entrySet())
		{
			if(seqName2node.containsKey(entry.getKey()))
			{
				assignmentNode2Clade.put(seqName2node.get(entry.getKey()), entry.getValue());
			}
			else
			{
				System.out.printf("WARNNING: given tree does not contain %s to assign %s%n", entry.getKey(), entry.getValue());
			}
		}
		
		ArrayList<FlexibleNode> removeNodes = new ArrayList<FlexibleNode>();
		for(HashMap.Entry<String, ArrayList<String>> entry : clade2samples.entrySet())
		{
			String clade = entry.getKey();
			ArrayList<String> samples = entry.getValue();
			ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();				
			for(int i=0; i<samples.size(); ++i)
			{
				if(seqName2node.containsKey(samples.get(i)))
				{
					FlexibleNode sampleNode = seqName2node.get(samples.get(i));
					ancestorList.clear(); get_parents2root(sampleNode, ancestorList);
					int j = 0;
					for(j=0; j<ancestorList.size(); ++j)
					{
						if(assignmentNode2Clade.containsKey(ancestorList.get(j)))
						{
							String predictClade = assignmentNode2Clade.get(ancestorList.get(j));
							if(!predictClade.equals(clade))
							{
								removeNodes.add(sampleNode);
								removedSampleAndClades.put(samples.get(i),clade);
							}
							break;
						}
					}
					if(j == ancestorList.size())
					{
						removeNodes.add(sampleNode);
						removedSampleAndClades.put(samples.get(i),clade);
					}
				}
				else
				{
					System.out.printf("WARNNING: given tree does not contain sample %s in clade/lineage %s%n", samples.get(i), clade);
				}
			}
		}
		//System.out.printf("INFO: remove %d from tree%n", removeNodes.size());
		MyFlexibleTree removedTree = removeTaxonReturnTree((FlexibleTree) tree, removeNodes);
		return removedTree;
	}
	
	private int getAmbiousNucleotidesNum(String seq)
	{
		int count = 0;
		for(int i=0; i<seq.length(); ++i)
		{
			char base = seq.charAt(i);
			if( _nucleotide_nomenclature.get((byte) base).size() > 1) count++;
		}
		return count;
	}
	
	private int getAmbiousNucleotidesNum(ConcurrentHashMap<Integer, Byte> seq)
	{
		int count = 0;
		for(HashMap.Entry<Integer, Byte> entry : seq.entrySet())
		{
			byte base = entry.getValue();
			if( _nucleotide_nomenclature.get(base).size() > 1) count++;
		}
		return count;
	}
	
	
	private ArrayList<String> sortSamplebyAmbiousCode(ArrayList<String> samples, String informat)
	{
		HashMap<String, Integer> seq2AmbiousCount = new HashMap<String, Integer>();
		for(int i=0; i<samples.size(); ++i)
		{
			int count = 0;
			if(informat.equals("vcf"))
			{
				ConcurrentHashMap<Integer, Byte> vcf_seq = getVariantSequenceByName(samples.get(i));
				count = getAmbiousNucleotidesNum(vcf_seq);
			}
			else
			{
				count = getAmbiousNucleotidesNum(getStringSequenceByName(samples.get(i)));
			}
			seq2AmbiousCount.put(samples.get(i), count);
		}
		
		seq2AmbiousCount = sortByValue(seq2AmbiousCount);
		ArrayList<String> sortedSamples = new ArrayList<String>(seq2AmbiousCount.keySet());
		return sortedSamples;
	}
	
	public Tree treeRefinement(Tree removedTree, HashMap<String, String> assignmentMap, HashMap<String, String> addSample2Clades, Integer index4Insertion, boolean printDisInfoOnScreen, String informat, String checkPoint_prefix) throws IOException, ImportException
	{
		return treeRefinement(removedTree, assignmentMap, addSample2Clades, index4Insertion, false, printDisInfoOnScreen, informat, checkPoint_prefix);
	}
	
	public Tree treeRefinement(Tree removedTree, HashMap<String, String> assignmentMap, HashMap<String, String> addSample2Clades, Integer index4Insertion, boolean stillInsert, boolean printDisInfoOnScreen, String informat, String checkPoint_prefix) throws IOException, ImportException
	{
		this.mytree = removedTree;
		setupHashtableOfnode2seqName();
		ArrayList<String> addSamples = new ArrayList<String>(addSample2Clades.keySet());
		ArrayList<String> sortedSamples = sortSamplebyAmbiousCode(addSamples, informat);
		
		HashMap<String, String> unAddSample2Clades = new HashMap<String, String>();
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();	
		int successfulPlacementNum = 0;

		for(int i=sortedSamples.size()-1; i>=0; --i)
		{
			String name = sortedSamples.get(i);
			System.out.println(i + ":" + name);
			try
			{
				Tree outtree = null;
				String qid = "q" + index4Insertion.intValue();
				String pid = "p" + index4Insertion.intValue();
				if(informat.contains("vcf"))
				{
					ConcurrentHashMap<Integer, Byte> query = getVariantSequenceByName(name);
					outtree = addQuerySequence(name, query, qid, pid, false, new double[3], "insertion_vcf", 0);
				}
				else
				{
					String query = getStringSequenceByName(name);
					outtree = addQuerySequence(name, query, qid, pid, false, new double[3], "insertion", 0);
				}
				
				
				index4Insertion++;
				
				//checkpoint to store the current outtree and remaining samples to be added
				if(checkPoint_prefix != "" && (i % 1000 == 999 || i == 0))
				{
					String outTreePath = OUTPUT_FOLDER  + "/" + checkPoint_prefix + "_addedTree.nwk";
					writeToTree(outtree, outTreePath, "refinement");
					String outFilePath = OUTPUT_FOLDER + "/" + checkPoint_prefix + "_remainingSeqName.txt";
					ArrayList<String> remainingSeqName = new ArrayList<String>();
					for(int j = 0; j < i; ++j) remainingSeqName.add(sortedSamples.get(j));
					writeSeqName(outFilePath, remainingSeqName);
					
					if(informat.contains("fasta"))
					{
						String outFASTAFilePath = OUTPUT_FOLDER + "/" + checkPoint_prefix + "_addedInternodeQ.fasta";
						ArrayList<String> desc = new ArrayList<String>();
						ArrayList<String> seq = new ArrayList<String>();
						for(HashMap.Entry<String, Integer> entry : seqIdxMap.entrySet()) {
							if(entry.getKey().charAt(0) == 'p')
							{
								desc.add(entry.getKey());
								seq.add(stringSequencesList.get(seqIdxMap.get(entry.getKey())));
							}
						}
						writeFASTA(desc, seq, outFASTAFilePath);
					}
				}
				
				
				if(stillInsert)  //add all samples to the tree
				{ 
					this.mytree = outtree;
					setupHashtableOfnode2seqName();
					continue;
				}
				
				//only keep correctly inserted sample in the tree
				HashMap<FlexibleNode, String> assignmentNode2Clade = new HashMap<FlexibleNode, String>();
				HashMap<String, FlexibleNode> seqName2node = setupHashtableOfseqName2node(outtree);
				for(HashMap.Entry<String,String> entry : assignmentMap.entrySet())
				{
					if(seqName2node.containsKey(entry.getKey()))
					{
						assignmentNode2Clade.put(seqName2node.get(entry.getKey()), entry.getValue());
					}
					else
					{
						System.out.printf("WARNNING: given tree does not contain %s to assign %s%n", entry.getKey(), entry.getValue());
					}
				}
				String clade = addSample2Clades.get(name);
				FlexibleNode sampleNode = seqName2node.get(name);
				ancestorList.clear(); get_parents2root(sampleNode, ancestorList);
				int j = 0;
				for(j=0; j<ancestorList.size(); ++j)
				{
					if(assignmentNode2Clade.containsKey(ancestorList.get(j)))
					{
						String predictClade = assignmentNode2Clade.get(ancestorList.get(j));
						if(!predictClade.equals(clade))
						{
							unAddSample2Clades.put(name, clade);
							if(printDisInfoOnScreen)
								System.out.printf("DEBUG: sample %s in %s clade/lineage has been assigned to a different %s clade/lineage%n", name, clade, predictClade);
						}
						else
						{
							this.mytree = outtree;
							setupHashtableOfnode2seqName();
							successfulPlacementNum++;
							//if(printDisInfoOnScreen)
							//	System.out.printf("sample %s in clade/lineage %s has been re-assigned successfully (%d)%n", name, clade, successfulPlacementNum);
						}
						break;
					}
				}
				if(j == ancestorList.size())
				{
					unAddSample2Clades.put(name, clade);
					if(printDisInfoOnScreen)
						System.out.printf("DEBUG: sample %s in clade/lineage %s has not been assigned to any clade/lineage%n", name, clade);
				}	
			}
			catch (Exception e) {
				System.out.printf("ERROR: sample %s can not be processed.%n", name);
				e.printStackTrace();
			}
		}
		if(!stillInsert)
		{
			System.out.printf("INFO: %d out of %d sampls have been assigned correctly%n", successfulPlacementNum, sortedSamples.size());
			addSample2Clades.clear();
			addSample2Clades.putAll(unAddSample2Clades);
		}
		else
		{
			System.out.printf("INFO: %d sampls have been added back to the tree%n", sortedSamples.size());
			addSample2Clades.clear();
		}
		//System.gc();
		return this.mytree;
	}
	
	public static HashMap<String, String> readAnnotations(String fn) {
		HashMap<String, String> nodename2clade = new HashMap<String, String>();

		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			String oneline = null;
			//// read header
			while ((oneline = br.readLine()) != null) {
				if (oneline.isEmpty())
					continue;
				String[] words = oneline.trim().split("\\s+");
				if (words.length == 2) {
					nodename2clade.put(words[0], words[1]);
				} else {
					System.out.println("WARNNING: get a wrong format clade sample " + oneline);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return nodename2clade;
	}
	
	public static HashMap<FlexibleNode, String> readAnnotations(String fn, HashMap<String, FlexibleNode> seqName2node) {
		//HashMap<String, String> nodename2clade = new HashMap<String, String>();
		HashMap<FlexibleNode, String> assignmentNode2Clade = new HashMap<FlexibleNode, String>();
		try {
			BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(fn)));
			String oneline = null;
			//// read header
			while ((oneline = br.readLine()) != null) {
				if (oneline.isEmpty())
					continue;
				String[] words = oneline.trim().split("\\s+");
				if (words.length == 2) {
					if(seqName2node.containsKey(words[0]))
						assignmentNode2Clade.put(seqName2node.get(words[0]), words[1]);
					else
						System.out.println("WARNNING: can not find assigned node " + words[0]);
				} else {
					System.out.println("WARNNING: get a wrong format clade sample " + oneline);
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
		return assignmentNode2Clade;
	}
	
	//keep all inconsistent samples in the tree
	public Tree doTreeRefinement(HashMap<String, ArrayList<String>> clade2samples, String informat, boolean printDisInfoOnScreen, HashMap<String, String> nodename2clade, String checkPoint_prefix)
	{
		if (nodename2clade == null) nodename2clade = treeAnnotation(clade2samples, 0, printDisInfoOnScreen); //tree annotation, F-Scoree = 0
		//Tree refTree = this.mytree.getCopy();
		MyFlexibleTree refTree = new MyFlexibleTree(this.mytree, true);
		copyAttributeFromOneNodeToAnother((FlexibleNode) this.mytree.getRoot(), (FlexibleNode) refTree.getRoot());

		HashMap<String, String> removedSampleAndClades = new HashMap<String, String> ();
		Tree removedTree = removeInconsistentSample(clade2samples, nodename2clade, this.mytree, removedSampleAndClades); //remove inconsistent taxa
		//for test
		String outTreePath = OUTPUT_FOLDER + removedTree.getExternalNodeCount() + "_removedTree.tree";
		writeToTree(removedTree, outTreePath, "refinement");
		String outUnKeepSamplesFile = OUTPUT_FOLDER + removedSampleAndClades.size() + "_unKeepSamples_" + removedTree.getExternalNodeCount() + "_tree.tsv";
		writeAssignClade(outUnKeepSamplesFile, removedSampleAndClades);
		
		Integer index4Insertion = 0;
		if(checkPoint_prefix != "")
		{
			String inTreePath = OUTPUT_FOLDER  + "/" + checkPoint_prefix + "_addedTree.nwk";
			Tree intree;
			try {
				NewickImporter intni = new NewickImporter(new FileReader(inTreePath));
				intree = intni.importTree(null);
				removedTree = intree;
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} catch (ImportException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			String inFilePath = OUTPUT_FOLDER + "/" + checkPoint_prefix + "_remainingSeqName.txt";
			ArrayList<String> innames = readSeqName(inFilePath);
		
			if(informat.contains("fasta"))
			{
				String inFASTAFilePath = OUTPUT_FOLDER + "/" + checkPoint_prefix + "_addedInternodeQ.fasta";
				readFastaAlignmentFile(inFASTAFilePath);
			}
			
			HashMap<String, String> current_removedSampleAndClades = new HashMap<String, String>();
			for(int i = 0; i < innames.size(); ++i)
			{
				current_removedSampleAndClades.put(innames.get(i), removedSampleAndClades.get(innames.get(i)));
			}
			removedSampleAndClades = current_removedSampleAndClades;
		}
		
		int current_max_used_index4Insertion = -1;
		for(HashMap.Entry<String, Integer> entry : seqIdxMap.entrySet()) {
			if(entry.getKey().charAt(0) == 'p')
			{
				int usedIdx = Integer.parseInt(entry.getKey().substring(1));
				if(usedIdx > current_max_used_index4Insertion) current_max_used_index4Insertion = usedIdx;
			}
		}
		index4Insertion = current_max_used_index4Insertion + 1;
		
		Tree outtree = null;
		try
		{
			outtree = treeRefinement(removedTree, nodename2clade, removedSampleAndClades, index4Insertion, true, printDisInfoOnScreen, informat, checkPoint_prefix); //tree refinement
		}
		catch (Exception e) {
			e.printStackTrace();
		}
		removedTree = outtree;

		if (outtree != null)
		{
			int totalTaxaInTree = refTree.getExternalNodeCount();
			int realTaxaInTree = outtree.getExternalNodeCount();
			if(totalTaxaInTree != realTaxaInTree)
			{
				System.out.println("ERROR: tree refinement with inconsistent taxa");
				System.exit(-1);
			}
		}
		return outtree;
	}
	
	public void getConsistentLabelTree(HashMap<String, ArrayList<String>> clade2samples, boolean printDisInfoOnScreen, HashMap<FlexibleNode, String> assignmentNode2Clade)
	{
		HashMap<String, String> nodename2clade = new HashMap<String, String>();
		for(HashMap.Entry<FlexibleNode, String> entry : assignmentNode2Clade.entrySet())
		{
			String selectNodeName = node2seqName.get(entry.getKey());
			nodename2clade.put(selectNodeName, entry.getValue());
		}
		
		if (nodename2clade == null) nodename2clade = treeAnnotation(clade2samples, 0, printDisInfoOnScreen); //tree annotation, F-Scoree = 0
		//Tree refTree = this.mytree.getCopy();
		MyFlexibleTree refTree = new MyFlexibleTree(this.mytree, true);
		copyAttributeFromOneNodeToAnother((FlexibleNode) this.mytree.getRoot(), (FlexibleNode) refTree.getRoot());

		HashMap<String, String> removedSampleAndClades = new HashMap<String, String> ();
		Tree removedTree = removeInconsistentSample(clade2samples, nodename2clade, this.mytree, removedSampleAndClades); //remove inconsistent taxa
		
		String outTreePath = OUTPUT_FOLDER + removedTree.getExternalNodeCount() + "_removedTree.nwk";
		writeToTree(removedTree, outTreePath, "refinement");
		String outUnKeepSamplesFile = OUTPUT_FOLDER + removedSampleAndClades.size() + "_unKeepSamples_" + removedTree.getExternalNodeCount() + "_tree.tsv";
		writeAssignClade(outUnKeepSamplesFile, removedSampleAndClades);
		
		System.out.println("INFO: remove inconsistent taxa and write this removed tree to file " + outTreePath);
	}
	
	public void getCollapsedTree(HashMap<FlexibleNode, String> assignmentNode2Clade)
	{
		//create a mapping to their ancestor
		HashMap<FlexibleNode, FlexibleNode> assignedNode2parent = new HashMap<FlexibleNode, FlexibleNode>();
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();
		for (HashMap.Entry<FlexibleNode, String> entry : assignmentNode2Clade.entrySet()) {
			//get the preceding annotated node
			FlexibleNode annotated_node = entry.getKey();
			ancestorList.clear(); get_parents2root(annotated_node, ancestorList);
			int j = 0;
			FlexibleNode annotated_node_precedor = ancestorList.get(ancestorList.size()-1); //default to be root
			for(j=1; j<ancestorList.size(); ++j) //the first item in ancestorList is annotated_node itself
			{
				if(assignmentNode2Clade.containsKey(ancestorList.get(j)))
				{
					annotated_node_precedor = ancestorList.get(j);
					break;
				}
			}
			assignedNode2parent.put(annotated_node, annotated_node_precedor);
		}
		HashMap<FlexibleNode, ArrayList<FlexibleNode>> parent2assignedNode = new HashMap<FlexibleNode, ArrayList<FlexibleNode>>();
		for (HashMap.Entry<FlexibleNode, FlexibleNode> entry : assignedNode2parent.entrySet()) {
			FlexibleNode child = entry.getKey();
			FlexibleNode parent = entry.getValue();
			if(!parent2assignedNode.containsKey(parent)) parent2assignedNode.put(parent, new ArrayList<FlexibleNode>());
			parent2assignedNode.get(parent).add(child);
		}

		FlexibleNode tree_root = (FlexibleNode) this.mytree.getRoot();
		for(int i = 0; i < parent2assignedNode.get(tree_root).size(); ++i) 
		{
			if(parent2assignedNode.get(tree_root).get(i) == tree_root)
			{
				parent2assignedNode.get(tree_root).remove(i);
				break;
			}
		}
		StringBuilder buffer = new StringBuilder();
		toMyNewick(parent2assignedNode, tree_root, 0, buffer, assignmentNode2Clade);
		String outTreePath = OUTPUT_FOLDER + assignmentNode2Clade.size() + "_collapsedTree.nwk";
		PrintStream out;
		try {
			out = new PrintStream(new FileOutputStream(new File(outTreePath)));
			out.println(buffer.toString());
			out.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		System.out.println("INFO: collapse the tree by lineages and write this collapsed tree to file " + outTreePath);
	}
	
	private void toMyNewick(HashMap<FlexibleNode, ArrayList<FlexibleNode>> parent2assignedNode, FlexibleNode current_node, double branch_len, StringBuilder buffer, HashMap<FlexibleNode, String> assignmentNode2Clade)
	{
		ArrayList<FlexibleNode> children = new ArrayList<FlexibleNode>();
		if(parent2assignedNode.containsKey(current_node)) children = parent2assignedNode.get(current_node);
		
		if(children.size() > 0) buffer.append('(');
		for(int i = 0; i < children.size(); ++i)
		{
			double child_parent_dist = 0;
			FlexibleNode my_parent = children.get(i);
			while(my_parent != current_node) 
			{
				child_parent_dist += my_parent.getLength();
				my_parent = my_parent.getParent();
			}
			toMyNewick(parent2assignedNode, children.get(i), child_parent_dist, buffer, assignmentNode2Clade);
			if(i < children.size() - 1) buffer.append(',');
		}
		if(children.size() > 0) buffer.append(')');
		String lineageName = "ROOT";
		if(assignmentNode2Clade.containsKey(current_node)) lineageName = assignmentNode2Clade.get(current_node);
		buffer.append(lineageName + "|" + node2seqName.get(current_node) + ":" + branch_len);
	}
	
	//only keep correct inserted samples in the tree that will do the TreeRefinement iteratively
	public Tree doTreeRefinement_iteration(HashMap<String, ArrayList<String>> clade2samples, String informat, boolean printDisInfoOnScreen, HashMap<String, String> nodename2clade, int iteration, String checkPoint_prefix)
	{
		if (nodename2clade == null) nodename2clade = treeAnnotation(clade2samples, 0, printDisInfoOnScreen); //tree annotation, F-Scoree = 0

		HashMap<String, String> removedSampleAndClades = new HashMap<String, String> ();
		MyFlexibleTree refTree = new MyFlexibleTree(this.mytree, true);
		copyAttributeFromOneNodeToAnother((FlexibleNode) this.mytree.getRoot(), (FlexibleNode) refTree.getRoot());
		
		HashMap<String, FlexibleNode> refSeqName2node = setupHashtableOfseqName2node(refTree); //only for test
		
		Tree removedTree = removeInconsistentSample(clade2samples, nodename2clade, this.mytree, removedSampleAndClades); //remove inconsistent taxa
		//for test
		String outTreePath = OUTPUT_FOLDER + removedTree.getExternalNodeCount() + "_removedTree.tree";
		writeToTree(removedTree, outTreePath, "refinement");
		String outUnKeepSamplesFile = OUTPUT_FOLDER + removedSampleAndClades.size() + "_unKeepSamples_" + removedTree.getExternalNodeCount() + "_tree.tsv";
		writeAssignClade(outUnKeepSamplesFile, removedSampleAndClades); ///write inconsistent sequences
		
		Integer index4Insertion = 0;
		Tree outtree = null;
		for(int i=0; i<iteration; ++i)
		{
			int addSamplesNum = removedSampleAndClades.size();
			try
			{
				//only keep correct inserted samples in the tree
				outtree = treeRefinement(removedTree, nodename2clade, removedSampleAndClades, index4Insertion, printDisInfoOnScreen, informat, checkPoint_prefix); //tree refinement
			}
			catch (Exception e) {
				e.printStackTrace();
			}
			removedTree = outtree;
			HashMap<String, ArrayList<String>> clade2samplesInTree = new HashMap<String, ArrayList<String>>();
			for(HashMap.Entry<String, ArrayList<String>> entry : clade2samples.entrySet())
			{
				ArrayList<String> sampleList = new ArrayList<String>();
				for(int j=0; j<entry.getValue().size(); ++j)
				{
					if(!removedSampleAndClades.containsKey(entry.getValue().get(j)))
					{
						sampleList.add(entry.getValue().get(j));
					}
				}
				if(sampleList.size() > 0)
				{
					clade2samplesInTree.put(entry.getKey(), sampleList);
				}
			}
		
			//only for test
			outUnKeepSamplesFile = OUTPUT_FOLDER + removedSampleAndClades.size() + "_unKeepSamples_" + outtree.getExternalNodeCount() + "_tree.tsv";
			writeAssignClade(outUnKeepSamplesFile, removedSampleAndClades);
			if (outtree != null)
			{
				int totalTaxaInTree = refTree.getExternalNodeCount() - removedSampleAndClades.size();
				int realTaxaInTree = outtree.getExternalNodeCount();
				if(totalTaxaInTree != realTaxaInTree)
				{
					System.out.println("ERROR: tree refinement with inconsistent taxa");
					System.exit(-1);
				}
				outTreePath = OUTPUT_FOLDER + "refined_" + realTaxaInTree + ".tree";
				writeToTree(outtree, outTreePath, "refinement");
				ArrayList<FlexibleNode> removedFromReftree = new ArrayList<FlexibleNode>();
				for(HashMap.Entry<String, String> entry : removedSampleAndClades.entrySet())
				{
					removedFromReftree.add(refSeqName2node.get(entry.getKey()));
				}
				Tree refRemovedTree = removeTaxonReturnTree((FlexibleTree) refTree, removedFromReftree);
				outTreePath = OUTPUT_FOLDER + "reference_" + realTaxaInTree + ".tree";
				writeToTree(refRemovedTree, outTreePath, "refinement");
				
				if(addSamplesNum == removedSampleAndClades.size()) //no improvement
				{
					///keep the inconsistent and unassigned samples in the tree
					try {
						outtree = treeRefinement(removedTree, nodename2clade, removedSampleAndClades, index4Insertion, true, printDisInfoOnScreen, informat, checkPoint_prefix);
					} catch (IOException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					} catch (ImportException e) {
						// TODO Auto-generated catch block
						e.printStackTrace();
					}
					break;
				}
			}
			else 
			{
				break;
			}
			
			if(addSamplesNum == removedSampleAndClades.size()) //no improvement
			{
				break;
			}
		}
		return outtree;
	}
	
	//to do in the future 
	public Tree doGraftSubtrees(HashMap<Integer, String> queryList, boolean printDisInfoOnScreen) throws InterruptedException
	{
		ArrayList<Integer> queryIdxs = new ArrayList<Integer>(queryList.keySet());
		queryIdxs.sort(Comparator.naturalOrder());
		placements = new String[queryIdxs.size()];
		Tree outtree = null;
		for (int i = 0; i < queryIdxs.size(); i++) {
			int idx = queryIdxs.get(i);
			String name = queryList.get(idx);
			String query = stringSequencesList.get(idx);
			String qid = "q" + (i + 1);
			String pid = "p" + (i + 1);
			outtree = addQuerySequence(name, query, qid, pid, printDisInfoOnScreen, new double[3], "placement", i);
		}
		
		int index = queryIdxs.size();
		for(HashMap.Entry<String, ArrayList<String>> entry : node2placementSeqName.entrySet())
		{
			String currentNodeBName = entry.getKey();
			ArrayList<String> queries =  entry.getValue();
			HashMap<String, FlexibleNode> seqName2node = setupHashtableOfseqName2node(this.mytree);
			if(!seqName2node.containsKey(currentNodeBName))
			{
				System.out.printf("ERROR: grafted node %s is not in reference tree", currentNodeBName);
				System.exit(-1);
			}
			FlexibleNode currentNodeB = seqName2node.get(currentNodeBName);
			
			if(queries.size() < 2)
			{
				//insert it back the reference tree directly
				String qname = queries.get(0);
				String nodeQseq = getStringSequenceByName(qname);
				String qid = "q" + (++index);
				String pid = "p" + (++index);
				outtree = addQuerySequence2Tree(currentNodeB, qname, nodeQseq, qid, pid, printDisInfoOnScreen);
				this.mytree = outtree;
				setupHashtableOfnode2seqName();
			}
			else
			{
				//call IQTREE2
				//exePath
				//step1. write the queries with A-B to file
				ArrayList<String> seqs = new ArrayList<String>();
				for(int i=0; i<queries.size(); ++i) 
					seqs.add(getStringSequenceByName(queries.get(i)));
				
				queries.add(currentNodeBName);
				seqs.add(getStringSequenceByName(currentNodeBName));
				String currentNodeAName = node2seqName.get(currentNodeB.getParent());
				
				queries.add(currentNodeAName);
				seqs.add(getStringSequenceByName(currentNodeAName));
				
				String output_folder = OUTPUT_FOLDER + currentNodeBName + "_attached_" + queries.size() + ".fasta";
				writeFASTA(queries, seqs, output_folder);
				
				//call iqtree2 to build up subtree
				String exePath = "";
				String params = "";
				ProcessBuilder pb = new ProcessBuilder(exePath, params);
				Process pro;
				try {
					pro = pb.start(); // Start the process.
					int exitVal = pro.waitFor();
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				} 
				
				//collect the subtree and reroot it by A
				
			}
		}
		
		return outtree;
	}
	
	public static Tree reRootbyGivenNode(Tree intree, String rerootName)
	{
		MyFlexibleTree mynewTree = new TIPars().new MyFlexibleTree(intree, true);
		copyAttributeFromOneNodeToAnother((FlexibleNode) intree.getRoot(), (FlexibleNode) mynewTree.getRoot());
		HashMap<String, FlexibleNode> seqName2node = setupHashtableOfseqName2node(mynewTree);
		FlexibleNode rerootNode = seqName2node.get(rerootName);
		
		if(rerootNode == null || rerootNode.isRoot()) return mynewTree;
			
		ArrayList<FlexibleNode> ancestorList = new ArrayList<FlexibleNode>();	
		get_parents2root(rerootNode, ancestorList);
		
		mynewTree.beginTreeEdit();
		double branchLen = 0;
		for(int i=ancestorList.size() - 1; i > 0; --i)
		{
			FlexibleNode set2Child = ancestorList.get(i);
			FlexibleNode set2Parent = ancestorList.get(i-1);
			if(i == ancestorList.size() - 1 && set2Child.getChildCount() == 2)
			{
				//deal with the original root with only two children
				set2Child = ancestorList.get(i).getChild(0);
				branchLen = ancestorList.get(i).getChild(1).getLength();
				if(set2Child.getNumber() == set2Parent.getNumber()) 
				{
					set2Child = ancestorList.get(i).getChild(1);
					branchLen = ancestorList.get(i).getChild(0).getLength();
				}
				branchLen += set2Child.getLength(); //merge the two sub-branches under the original root 
				set2Child.setLength(branchLen);
			}
			else
			{
				branchLen = set2Parent.getLength();
				set2Child.setLength(branchLen);
				set2Child.removeChild(set2Parent);
			}
			set2Child.setParent(set2Parent);
			set2Parent.addChild(set2Child);
		}
		mynewTree.setRoot(rerootNode);
		mynewTree.endTreeEdit();
		mynewTree.toAdoptNodes((FlexibleNode) mynewTree.getRoot());
		
		return mynewTree;	
	}

	public static class JSONLink {
		String source;
		String target;
		double distance;

		public JSONLink(String source, String target, double distance) {
				this.source = source;
				this.target = target;
				this.distance = distance;
		}
	}

	public static class JSONNode {
		String id;
		int numberOfSamples;
		boolean isRoot = false;
		String cladeRoot;
		String annotation;
		double distanceToRoot;
		int numberOfTruePositive;
		int numberOfFalsePositive;
		int numberOfFalseNegative;
		int numberOfUnlabeledSamples;
		double f1score;

		public JSONNode(
			String id,
			int numberOfSamples,
			boolean isRoot,
			String cladeRoot,
			String annotation,
			double distanceToRoot,
			int numberOfTruePositive,
			int numberOfFalsePositive,
			int numberOfFalseNegative,
			int numberOfUnlabeledSamples,
			double f1score
		) {
			this.id = id;
			this.numberOfSamples = numberOfSamples;
			this.isRoot = isRoot;
			this.cladeRoot = cladeRoot;
			this.annotation = annotation;
			this.distanceToRoot = distanceToRoot;
			this.numberOfTruePositive = numberOfTruePositive;
			this.numberOfFalsePositive = numberOfFalsePositive;
			this.numberOfFalseNegative = numberOfFalseNegative;
			this.numberOfUnlabeledSamples = numberOfUnlabeledSamples;
			this.f1score = f1score;
		}
	
		public JSONNode(
			String id,
			int numberOfSamples,
			boolean isRoot
		) {
			this.id = id;
			this.numberOfSamples = numberOfSamples;
			this.isRoot = isRoot;
		}
	}

	public static class JSONResult {
		ArrayList<JSONLink> links;
		ArrayList<JSONNode> nodes;
		public JSONResult() {
			this.links = new ArrayList<>();
			this.nodes = new ArrayList<>();
		}
	}

	public static class CladeDetailLineParser {
		private HashMap<String, String> sample2Clade;
		private HashMap<String, Integer> clade2NumberOfSamples;
		public JSONLink resultLink;
		public JSONNode resultNode;
		public String rootNodeName = null; // for root without annotation ONLY
		public int totalSamplesCount = 0;
	
		public CladeDetailLineParser(
			HashMap<String, ArrayList<String>> clade2Samples
		) {
			setHashMaps(clade2Samples);
		}

		private void setHashMaps(
			HashMap<String, ArrayList<String>> clade2Samples
		) {
			HashMap<String, String> sample2Clade = new HashMap<>();
			HashMap<String, Integer> clade2NumberOfSamples = new HashMap(clade2Samples.size());
			for (String clade : clade2Samples.keySet()) {
				ArrayList<String> samples = clade2Samples.get(clade);
				for (String sample : samples) {
					sample2Clade.put(sample, clade);
				}
				int sampleSize = samples.size();
				clade2NumberOfSamples.put(clade, sampleSize);
				this.totalSamplesCount += sampleSize;
			}
			this.clade2NumberOfSamples = clade2NumberOfSamples;
			this.sample2Clade = sample2Clade;
		}
		
		public void parseLine(String line) throws Exception {
			String[] words = line.split("\t");
			if (words.length != 6) {
				throw new Exception("words length not equals 6 : " + line);
			}
			final boolean isRoot = "NA".equals(words[1]);
			final boolean isAnnotated = !"null".equals(words[3]);
			if (isRoot && !isAnnotated) {
				this.rootNodeName = words[0];
				return;
			}
			final double distanceToRoot = Double.parseDouble(words[2]);;
			final double f1score = Double.parseDouble(words[4]);
			final String annotation = words[3];
			String[] samplesInClade = words[5].split(",");
			int samplesCount = samplesInClade.length,
				truePositiveCount = 0, 
				falsePositiveCount = 0, 
				unlabeledCount = 0;
			for (String sample : samplesInClade) {
				if(sample.isEmpty() || sample.isBlank()) continue;
				String clade = sample2Clade.get(sample);
				if (clade == null) {
					unlabeledCount++;
				} else if (clade.equals(annotation)) {
					truePositiveCount++;
				} else {
					falsePositiveCount++;
				}
			}
			if (truePositiveCount + falsePositiveCount + unlabeledCount != samplesCount) {
				throw new Exception(String.format(
					"true positive count(%d) + false positive count(%d) + unlabeled count(%d) != smaples Count(%d)",
					truePositiveCount,
					falsePositiveCount,
					unlabeledCount,
					samplesCount	
				));
			}
			final int numberOfFalseNegative = clade2NumberOfSamples.get(annotation) - truePositiveCount;
			resultNode = new JSONNode(
				words[0],
				samplesCount,
				isRoot,
				words[0],
				words[3],
				distanceToRoot,
				truePositiveCount,
				falsePositiveCount,
				numberOfFalseNegative,
				unlabeledCount,
				f1score
			);
			if (!isRoot) { // is not root
				resultLink = new JSONLink(
					words[1],
					words[0],
					distanceToRoot
				);
			}
		}
	}

	public static void generateGraphJSON(
		String jsonOutputFilename,
		String annotationTooltipName,
		String cladesDetailsFilename,
		HashMap<String, ArrayList<String>> clade2samples
	) {
		Path cladesDetailsFilePath = Path.of(cladesDetailsFilename);
		CladeDetailLineParser parser = new CladeDetailLineParser(clade2samples);
		final int numberOfClades = clade2samples.size();
		JSONResult jsonResult = new JSONResult();
		try (Stream<String> lines = Files.lines(cladesDetailsFilePath)) {
			int lineCount = 0;
			for (String line : (Iterable<String>) lines::iterator) {
				if (lineCount > 0) { // skip column name line
					try {
						parser.parseLine(line);
						jsonResult.links.add(parser.resultLink);
						jsonResult.nodes.add(parser.resultNode);	
					} catch(Exception e) {
						e.printStackTrace();
					}
				}
				lineCount++;
			}
			if (parser.rootNodeName != null) { // root without annotation
				JSONNode rootClade = new JSONNode(
					parser.rootNodeName,
					parser.totalSamplesCount,
					true
				);
				jsonResult.nodes.add(rootClade);
			}
			//System.out.println("node size: " + jsonResult.nodes.size());
			//System.out.println("link size: " + jsonResult.links.size());
			Path jsonOutputFilePath = Path.of(jsonOutputFilename);
			GsonBuilder gsonBuilder = new GsonBuilder();
        	gsonBuilder.setPrettyPrinting();
			Gson gson = gsonBuilder.create();
			String json = gson.toJson(jsonResult);
			String fileData = "const data = " + json + ";\nconst annotationTooltipName = '" + annotationTooltipName + "';";
			Files.writeString(jsonOutputFilePath, fileData, StandardCharsets.UTF_8);
			//System.out.println("JS graph data files is  written successfully, path: " + jsonOutputFilename);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public static void copyStaticWebFiles(String output_dir) {
		final String currentJarPath = TIPars.class.getProtectionDomain().getCodeSource().getLocation().getPath();
		final File jarFile = new File(currentJarPath);
		final String jarDir = jarFile.getParentFile().getAbsolutePath();
		final String refinedOutputPath = output_dir.isEmpty() ? "." : output_dir;
		try {
			copyFiles(jarDir + "/visual/graph.html", refinedOutputPath + "/graph.html");
			copyFiles(jarDir + "/visual/graph.css", refinedOutputPath + "/graph.css");
			copyFiles(jarDir + "/visual/graph.js", refinedOutputPath + "/graph.js");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static void copyFiles(String originalPath, String targetPath) {
		try {
			FileInputStream fis = new FileInputStream(originalPath);
			FileOutputStream fos = new FileOutputStream(new File(targetPath));
            byte[] buffer = new byte[1024];
            int bytesRead;
            while ((bytesRead = fis.read(buffer)) != -1) {
                fos.write(buffer, 0, bytesRead);
            }
        } catch (IOException e) {
            System.err.println("Failed to copy files: " + e.getMessage());
        }
	}
	
	public static String changeFileExtension(String fileName, String newExtension) {
		int dotIndex = fileName.lastIndexOf('.');
        
        if (dotIndex != -1) {
            String baseName = fileName.substring(0, dotIndex);
            return baseName + "." + newExtension;
        }
        
        return fileName + "." + newExtension;
	}

	
	private static HashMap<FlexibleNode, HashMap<String, HashSet<String>>> stratifyTree(Tree mytree, HashMap<FlexibleNode, ScoreAndStringList> id2PajekVectix, int exploreTreeNodeLimit, int smallBubbleLimit, int smallClusterLimit, HashMap<String, FlexibleNode> seqName2node) 
	{
		HashMap<FlexibleNode, HashMap<String, HashSet<String>>> stratification4eachVectix = new HashMap<FlexibleNode, HashMap<String, HashSet<String>>>();
		HashMap<FlexibleNode, ScoreAndStringList> small_id2PajekVectix = new HashMap<FlexibleNode, ScoreAndStringList>();
		HashMap<FlexibleNode, Integer> internalNodeSeqName2Idx = new HashMap<FlexibleNode, Integer>();
		HashMap<FlexibleNode, Integer> externalNodeSeqName2Idx = new HashMap<FlexibleNode, Integer>();
		for (int i = 0; i < mytree.getInternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) mytree.getInternalNode(i);
			internalNodeSeqName2Idx.put(n, i);
		}
		for (int i = 0; i < mytree.getExternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) mytree.getExternalNode(i);
			externalNodeSeqName2Idx.put(n, i);
		}
		
		//entry.getValue().nodeList.size() >= smallBubbleLimit
		for(Entry<FlexibleNode, ScoreAndStringList> entry : id2PajekVectix.entrySet())
		{
			if(entry.getValue().SampleList.size() > smallClusterLimit)
			{
				HashMap<String, HashSet<String>> stratification = stratifyTreeByNode(mytree, entry, exploreTreeNodeLimit, smallBubbleLimit, seqName2node, internalNodeSeqName2Idx, externalNodeSeqName2Idx);
				if(stratification != null) stratification4eachVectix.put(entry.getKey(), stratification);
			}
			else
			{
				small_id2PajekVectix.put(entry.getKey(), entry.getValue());
			}
		}
		
		//merge small_id2PajekVectix (un-stratify) to stratification4eachVectix (already stratify)
		//based on the minimum distance to the already setted bubbles which are their ancestors
		for(Entry<FlexibleNode, ScoreAndStringList> entry : small_id2PajekVectix.entrySet())
		{
			String bubbleId = node2seqName.get(entry.getKey());
			FlexibleNode bubble_node = entry.getKey();
			FlexibleNode bubble_node_setted_select_vertixId = null;
			String bubble_node_setted_select_bubbleId = null;
			double min_distance = 10000;
			
			ArrayList<String> ancestors = getAncestors(bubble_node, (FlexibleNode)mytree.getRoot());
			ancestors.add(node2seqName.get((FlexibleNode)mytree.getRoot()));
			HashSet<String> ancestors_set = new HashSet<String>(ancestors);
			
			for(Entry<FlexibleNode, HashMap<String, HashSet<String>>> entry1 : stratification4eachVectix.entrySet())
			{
				if(entry1.getValue() == null) continue;
				for(Entry<String, HashSet<String>>  entry2 : entry1.getValue().entrySet())
				{
					String bubbleId_setted = entry2.getKey();
					if(!ancestors_set.contains(bubbleId_setted)) continue;
					
					FlexibleNode bubble_node_setted = seqName2node.get(bubbleId_setted);
					FlexibleNodeBranch nb = returnCommonAncestor(bubble_node, bubble_node_setted, mytree);
					if(nb.b < min_distance)
					{
						min_distance = nb.b;
						bubble_node_setted_select_vertixId = entry1.getKey();
						bubble_node_setted_select_bubbleId = bubbleId_setted;
					}
				}
			}
			if(bubble_node_setted_select_vertixId != null) 
			{
				stratification4eachVectix.get(bubble_node_setted_select_vertixId).get(bubble_node_setted_select_bubbleId).addAll(entry.getValue().SampleList);
				//System.out.println("processing " + bubbleId + " with #taxa " + entry.getValue().SampleList.size() + ", merged into " + bubble_node_setted_select_bubbleId);
			}
			else
			{
				//can not merge entry.
				//System.out.println("can not merge " + node2seqName.get(entry.getKey()));
				HashMap<String, HashSet<String>> temp = new HashMap<String, HashSet<String>>();
				temp.put(node2seqName.get(entry.getKey()), new HashSet<String>(entry.getValue().SampleList));
				stratification4eachVectix.put(entry.getKey(), temp);
				//System.out.println("processing " + bubbleId + " with #taxa " + entry.getValue().SampleList.size());
			}
		}

		return stratification4eachVectix;
	}
	
	private static HashMap<String, HashSet<String>> stratifyTreeByNode(Tree mytree, Entry<FlexibleNode, ScoreAndStringList> vertix, int exploreTreeNodeLimit, int smallBubbleLimit, HashMap<String, FlexibleNode> seqName2node, HashMap<FlexibleNode, Integer> internalNodeSeqName2Idx, HashMap<FlexibleNode, Integer> externalNodeSeqName2Idx)
	{
		FlexibleNode bubbleId_node = vertix.getKey();
		String bubbleId = node2seqName.get(bubbleId_node);
		
		if(vertix.getValue().SampleList.size() == 1)
		{
			HashMap<String, HashSet<String>> one_stratification = new HashMap<String, HashSet<String>>();
			HashSet<String> temp = new HashSet<String>();
			temp.add(bubbleId);
			one_stratification.put(bubbleId, temp);
			return one_stratification;
		}
			
		MyFlexibleTree reroottree = new TIPars().new MyFlexibleTree(mytree, true);
		copyAttributeFromOneNodeToAnother((FlexibleNode) mytree.getRoot(), (FlexibleNode) reroottree.getRoot());
		FlexibleNode newroot = null;
		
		if(!bubbleId_node.isExternal()) 
		{
			int internalIdx = internalNodeSeqName2Idx.get(bubbleId_node);
			newroot = (FlexibleNode) reroottree.getInternalNode(internalIdx);
		}
		else
		{
			int externalIdx = externalNodeSeqName2Idx.get(bubbleId_node);
			newroot = (FlexibleNode) reroottree.getExternalNode(externalIdx);
		}

		reroottree.beginTreeEdit();
		if(!newroot.isRoot()) reroottree.setRoot(newroot);
		reroottree.endTreeEdit();
		reroottree.toAdoptNodes(newroot);
		
		HashMap<FlexibleNode, String> myNode2seqName = setupNode2seqName(reroottree);
		
		//regenerate a subtree contains all vertix.nodelist 
		ArrayList<FlexibleNode> removeNodes = new ArrayList<FlexibleNode>();
		for(int i = 0; i < reroottree.getExternalNodeCount(); ++i)
		{
			FlexibleNode node = (FlexibleNode) reroottree.getExternalNode(i);
			
			if(!vertix.getValue().SampleList.contains(myNode2seqName.get(node))) 
			{
				removeNodes.add(node);
			}
		}
				
		Tree prundedtree = removeTaxonReturnTree((FlexibleTree) reroottree, removeNodes);
		
		for (int i = 0; i < prundedtree.getExternalNodeCount(); i++) {
			FlexibleNode n = (FlexibleNode) prundedtree.getExternalNode(i);
			String sequenceName = n.getTaxon().getId();
			//if(!vertix.getValue().SampleList.contains(sequenceName)) System.out.println(sequenceName + " not found");
		}
		
		//System.out.println(reroottree.getExternalNodeCount() + "-" + removeNodes.size() + " = " + prundedtree.getExternalNodeCount());
		//System.out.println("processing " + bubbleId + " with #taxa " + vertix.getValue().SampleList.size());
		
		myNode2seqName = setupNode2seqName(prundedtree);
		
		HashMap<String, HashSet<String>> stratification = stratifyTree(prundedtree, exploreTreeNodeLimit, myNode2seqName);
		
		stratification = mergeSmallBubble(stratification, smallBubbleLimit);

		return stratification;
	}
	
	//get all parent nodes, excluding child and parent
    private static ArrayList<String> getAncestors(FlexibleNode child, FlexibleNode parent)
    {
    	ArrayList<String> ancestors = new ArrayList<String>();
    	FlexibleNode node = child.getParent();
    	while(node != null && node != parent)
    	{
    		ancestors.add(node2seqName.get(node));
    		node = node.getParent();
    	}
    	if(node == null)
    	{
    		//System.out.println(node2seqName.get(parent) + " is not parent of " + node2seqName.get(child));
    		return null;
    	}
    	return ancestors;
    }
	
  //get all parent nodes to root
    private static ArrayList<String> getAncestors(FlexibleNode child)
    {
    	ArrayList<String> ancestors = new ArrayList<String>();
    	FlexibleNode node = child.getParent();
    	while(node != null)
    	{
    		ancestors.add(node2seqName.get(node));
    		if(node.isRoot()) break;
    		else node = node.getParent();
    	}
  
    	return ancestors;
    }
	
    private static HashMap<String, HashSet<String>> stratifyTree(Tree mytree, int numLimit, HashMap<FlexibleNode, String> mynode2seqName) {
		HashMap<String, HashSet<String>> stratification = new HashMap<String, HashSet<String>>();
		FlexibleNode node = (FlexibleNode) mytree.getRoot();
		Queue<FlexibleNode> queue = new LinkedList<FlexibleNode>();
		queue.add(node);
		while (!queue.isEmpty()) {
			node = queue.poll();

			ExploreNodes explorenodes = bfsSearch(mytree, node, numLimit, mynode2seqName);
			stratification.put(mynode2seqName.get(node), explorenodes.exploreList);
			for (int i = 0; i < explorenodes.bubbleList.size(); ++i) {
				queue.add(explorenodes.bubbleList.get(i));
			}
		}
		return stratification;
	}
    
    static Comparator<FlexibleNode> comparatorChildCount = new Comparator<FlexibleNode>() {
		public int compare(FlexibleNode p1, FlexibleNode p2) {
			if (!p1.hasChildren())
				return -1;
			else if (!p2.hasChildren())
				return 1;

			if (p1.getChildCount() > p2.getChildCount())
				return 1;
			else if (p1.getChildCount() < p2.getChildCount())
				return -1;
			else
				return 0;
		}
	};
	
	class ExploreNodes {
		ExploreNodes(HashSet<String> _exploreList, ArrayList<FlexibleNode> _bubbleList) {
			exploreList = _exploreList;
			bubbleList = _bubbleList;
		}

		HashSet<String> exploreList;
		ArrayList<FlexibleNode> bubbleList;
	}
	
	private static ExploreNodes bfsSearch(Tree mytree, FlexibleNode myroot, int numlimit, HashMap<FlexibleNode, String> mynode2seqName) {

		HashSet<String> exploreList = new HashSet<String>();
		ArrayList<FlexibleNode> bubbleList = new ArrayList<FlexibleNode>();
		ArrayList<FlexibleNode> childNodeList = new ArrayList<FlexibleNode>();
		Queue<FlexibleNode> queue = new LinkedList<FlexibleNode>();
		queue.add(myroot);
		FlexibleNode node, childnode;
		int count = 1, nonzeroBrCount = 0; // set count = 1 because myroot is in.
		double branchLength;
		while (!queue.isEmpty()) {
			node = queue.peek();
			if (node.getChildCount() + count > numlimit && count > 1) // in case of 1 + myroot.childCount > numlimit
			{
				// stop BFS and set bubbleList
				while (!queue.isEmpty()) {
					node = queue.poll();
					if (node.getChildCount() > 0)
						bubbleList.add(node);
				}
				break;
			}
			nonzeroBrCount = 0;
			childNodeList.clear();
			for (int i = 0; i < node.getChildCount(); ++i) {
				if(node.getChild(i) != null) childNodeList.add(node.getChild(i));
			}
			// sort the childNodeList
			childNodeList.sort(comparatorChildCount);
			for (int i = 0; i < childNodeList.size(); ++i) {
				childnode = childNodeList.get(i);
				branchLength = childnode.getLength();
				if (!exploreList.contains(mynode2seqName.get(childnode))) {
					if (childnode.hasChildren())
						queue.add(childnode);
					exploreList.add(mynode2seqName.get(childnode));
					if (branchLength > 0.0000000001)
						nonzeroBrCount++; /// new add on 20220314, do not count zero branch
				}
			}
			count += nonzeroBrCount;
			queue.poll();
		}

		ExploreNodes explorenodes = new TIPars().new ExploreNodes(exploreList, bubbleList);
		return explorenodes;
	}

	public static HashMap<String, HashSet<String>> mergeSmallBubble(HashMap<String, HashSet<String>> old_stratification, int smallBubbleLimit)
	{
		if(old_stratification.size() == 1) return old_stratification;
		
		HashMap<String, HashSet<String>> new_stratification = new HashMap<String, HashSet<String>>();
		ArrayList<String> smallBubbleList = new ArrayList<String>();
		for(HashMap.Entry<String, HashSet<String>> entry : old_stratification.entrySet()) {
			if(entry.getValue().size() > smallBubbleLimit) new_stratification.put(entry.getKey(),entry.getValue());
			else smallBubbleList.add(entry.getKey());
		}
		
		for(int i = 0; i < smallBubbleList.size(); ++i)
		{
			String bubbleName = smallBubbleList.get(i);
			boolean isfound = false;
			for(HashMap.Entry<String, HashSet<String>> entry : new_stratification.entrySet())
			{
				if(entry.getValue().contains(bubbleName))
				{
					entry.getValue().addAll(old_stratification.get(bubbleName));
					isfound = true;
				}
			}
			if(!isfound) 
			{
				for(int j = i + 1; j < smallBubbleList.size(); ++j) //is in the other smallBubbleList
				{
					String othersmallbubbleName = smallBubbleList.get(j);
					if(old_stratification.get(othersmallbubbleName).contains(bubbleName))
					{
						old_stratification.get(othersmallbubbleName).addAll(old_stratification.get(bubbleName));
						isfound = true;
					}
				}
			}
			if(!isfound)
			{
				System.out.println("can not find small bubble " + bubbleName);
			}
		}
		return new_stratification;
	}
	
	/*
	 * three types of bubbles
	 * 1: the first (highest) level, output by TIPars2 annotation, # > smallClusterLimit
	 * 2: the second level, stratify for those large bubbles in level1 (# > exploreTreeNodeLimit)
	 */
	public static void writeAllStratification2TSV(HashMap<FlexibleNode, HashMap<String, HashSet<String>>> allstratification,  HashMap<FlexibleNode, ScoreAndStringList> id2PajekVectix, HashMap<FlexibleNode, String> assignedNode2Clade, HashMap<String, FlexibleNode> seqName2node, String outfile) 
	{
		try {
			StringBuilder buffer = new StringBuilder();
			buffer.append("bubble_type\tannotated_node\tannotated_node_precedor\tdist_to_precedor\tparent_node\tpangolineage\tnum_nodes\tnodes\n");
			HashSet<FlexibleNode> allHighestBubbles = new HashSet<FlexibleNode>(allstratification.keySet());
			ArrayList<FlexibleNode> allAncestors = new ArrayList<FlexibleNode>();
			for(Entry<FlexibleNode, HashMap<String, HashSet<String>>> allentry : allstratification.entrySet())
			{
				FlexibleNode pajekVectixIdx = allentry.getKey();
				String bubbleId = node2seqName.get(pajekVectixIdx);
				allAncestors.clear(); get_parents2root(pajekVectixIdx, allAncestors);
				FlexibleNode preAncestor = allAncestors.get(allAncestors.size()-1); //default to be root
				//distance to root
				double dist_to_preAncestor = 0;
				for(int i=1; i < allAncestors.size(); ++i) //exclude itself
				{
					if(allHighestBubbles.contains(allAncestors.get(i)))
					{
						preAncestor = allAncestors.get(i);
						break;
					}
					dist_to_preAncestor += allAncestors.get(i).getLength();
				}
				String top_bubbleId = bubbleId;
				
				//the parent of root is set to be null
				buffer.append(1 + "\t" + top_bubbleId + "\t" + (preAncestor != pajekVectixIdx? node2seqName.get(preAncestor) : "null") + "\t" + dist_to_preAncestor);
				
                if (!pajekVectixIdx.isRoot())	buffer.append("\t" + node2seqName.get(pajekVectixIdx.getParent())); //parent_node
                else buffer.append("\tnull");
                
                String ancestral_state = assignedNode2Clade.get(pajekVectixIdx);
                buffer.append("\t" + ancestral_state); //pangolineage
                
				HashMap<String, HashSet<String>> stratification = allentry.getValue();
				
				HashSet<String> top_bubbleId_nodelist = getBubbleSubTreeAllNodesInFullTree(stratification.get(top_bubbleId), top_bubbleId, seqName2node);
				buffer.append("\t" + stratification.get(top_bubbleId).size() + "\t");
					
				int count = 0;
				for (String nodeName : top_bubbleId_nodelist)
				{
					count++;
					if(count < top_bubbleId_nodelist.size())
						buffer.append(nodeName + ",");
					else
						buffer.append(nodeName + "\n");
				}
					
				if(stratification.size() > 1)
				{						
					for (Entry<String, HashSet<String>> entry : stratification.entrySet()) {
						if(entry.getKey().compareTo(top_bubbleId) != 0)
						{
							FlexibleNode annotated_node = seqName2node.get(entry.getKey());
							allAncestors.clear(); get_parents2root(annotated_node, allAncestors);
			
							dist_to_preAncestor = 0;
							FlexibleNode annotated_node_preAncestor = allAncestors.get(allAncestors.size()-1); //default to be root
							for(int i=1; i < allAncestors.size(); ++i) //exclude itself
							{
								if(stratification.containsKey(node2seqName.get(allAncestors.get(i))))
								{
									annotated_node_preAncestor = allAncestors.get(i);
									break;
								}
								dist_to_preAncestor += allAncestors.get(i).getLength();
							}
							buffer.append(2 + "\t" + entry.getKey() + "\t" + ((annotated_node_preAncestor != annotated_node)? node2seqName.get(annotated_node_preAncestor) : "null") + "\t" + dist_to_preAncestor);
							if (!annotated_node.isRoot())	buffer.append("\t" + node2seqName.get(annotated_node.getParent())); //parent_node
			                else buffer.append("\tnull");
							
							buffer.append("\t" + ancestral_state); //pangolineage
							
							HashSet<String> nodelist = getBubbleSubTreeAllNodesInFullTree(stratification.get(entry.getKey()), entry.getKey(), seqName2node);
							buffer.append("\t" + entry.getValue().size() + "\t");
							count = 0;
							for (String nodeName : nodelist)
							{
								count++;
								if(count < nodelist.size())
									buffer.append(nodeName + ",");
								else
									buffer.append(nodeName + "\n");
							}
						}
					}
				}
			}
			PrintStream out = new PrintStream(new FileOutputStream(new File(outfile)));
			out.println(buffer.toString());
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static HashSet<String> getBubbleSubTreeAllNodesInFullTree(HashSet<String> subtreeNodes, String subtree_root, HashMap<String, FlexibleNode> seqName2node)
	{
		HashSet<String> allSubtreeNodes = (HashSet<String>)subtreeNodes.clone();
		FlexibleNode root_node = seqName2node.get(subtree_root);
		for(String nodename : subtreeNodes)
		{
			FlexibleNode child_node = seqName2node.get(nodename);
			ArrayList<String> ancestors = getAncestors(child_node, root_node);
			if(ancestors != null)
			{
				for(int i=0; i<ancestors.size(); ++i)
				{
					if(!allSubtreeNodes.contains(ancestors.get(i)))
						allSubtreeNodes.add(ancestors.get(i));
				}
			}
		}
		return allSubtreeNodes;
	}
	

	
	public static void main(String[] args) {
		
		Runtime run = Runtime.getRuntime();
		long startTime = System.currentTimeMillis();
		
		String insfn = "";
		String intfn = "";
		String inafn = "";
		String inqfn = "";
		String informat = "fasta";
		String outfn = "";
		String otype = "insertion";
		Boolean aa_flag = false;
		double fscore_min = 0;
		String checkpoint_prefix = ""; //only available when using fasta file format for insertion of a set of sequences
		String annotation_assignment = "";
		//tree BFS
		int exploreTreeNodeLimit = 2000;
		int smallBubbleLimit = 100;
		int smallClusterLimit = 10;
		
		String nidname = "label";
		String attname = "GenName";
		boolean printDisInfoOnScreen = true;

		_nucleotide_nomenclature = generate_nucleotide_nomenclature();  // IUPAC nucleotide codes
		_nucleotide_nomenclature_scoreTable = generate_nucleotide_nomenclature_scoreTable(); // substitution table																						
		_nucleotide_nomenclature_map2char = generate_nucleotide_nomenclature_characterTable();
		
		_aminoacid_scoreTable = generate_aminoacid_scoreTable(); //Amino acid substitution table using blosum62 matrix

		try {
//			intfn = "/home/d24h_prog2/yytao/96020_seq/100K_Short_Br_removed.tree"; // args[0];
//			//intfn = "/home/d24h_prog2/yytao/96020_seq/annotation/AnotationAfterRemove/81075_removedTree.tree"; // args[0];
//			insfn = "/home/d24h_prog2/yytao/96020_seq/combined_ancestral_seq.fasta.taxa"; // args[1];
//			inafn = "/home/d24h_prog2/yytao/96020_seq/combined_ancestral_seq.fasta.ancseq"; // args[2];
//			//inqfn = "/home/d24h_prog2/yytao/96020_seq/annotation/AnotationAfterRemove/81075_pangolin.tsv"; // args[3];
//			inqfn = "/home/d24h_prog2/yytao/96020_seq/annotation/96020_pangolin.tsv"; // args[3];
//			outfn = "/home/d24h_prog2/yytao/96020_seq/annotation/test/refined_96020.tree"; // args[6];

//			intfn = "/home/d24h_prog2/yytao/gisaid_20210906/for_TIPars/global_with_innode_root.tree"; // args[0];
//			insfn = "/home/d24h_prog2/yytao/gisaid_20210906/for_TIPars/combined_all.taxa"; // args[1];
//			inafn = "/home/d24h_prog2/yytao/gisaid_20210906/for_TIPars/combined_all.ancseq"; // args[2];
//			inqfn = "/home/d24h_prog2/yytao/gisaid_20210906/pangolin/656264_pangolin";
//			//inqfn = "/home/d24h_prog2/yytao/gisaid_20210906/pangolin/648325_filter50_pangolin"; // args[3];
//			outfn = "/home/d24h_prog2/yytao/gisaid_20210906/pangolin/tipars2/660k_refined.tree"; // args[6];
			
			
			otype = args[0]; //"insertion","placement","refinement","refinement_from_annotation","annotation","annotation_details"
			intfn = args[1]; //treefile

			if(otype.equals("insertion") || otype.equals("placement"))
			{
	        	insfn = args[2]; //taxa
	            inafn = args[3]; //ances
	            inqfn = args[4]; //query
	            informat = args[5];   //format
	            isMultiplePlacements = Boolean.parseBoolean(args[6]); //is multiple placements or not
			}

			if(otype.equals("annotation"))
			{
				inqfn = args[2]; //clade2samples, the annotation table
				fscore_min = Double.parseDouble(args[3]);
			}
			
			
			if(otype.equals("annotation_details"))
			{
				inqfn = args[2]; //clade2samples, the annotation table
				fscore_min = 0;
				annotation_assignment = args[3]; //annotation_assignment file
			}
			
			if(otype.equals("refinement") || otype.equals("refinement_from_annotation"))
			{
				insfn = args[2]; //taxa
	            inafn = args[3]; //ances
	            inqfn = args[4]; //clade2samples, the annotation table
	            informat = args[5]; //format
	            
	            isMultiplePlacements = Boolean.FALSE;
				fscore_min = 0;
				
				if(otype.equals("refinement_from_annotation")) annotation_assignment = args[6]; 
			}
			
			if(otype.equals("tree_BFS"))
			{
				inqfn = args[2]; //clade2samples, the annotation table
				exploreTreeNodeLimit = Integer.parseInt(args[3]);
				smallBubbleLimit = Integer.parseInt(args[4]);
				smallClusterLimit = Integer.parseInt(args[5]);
			}
				 
			
			if (otype.equals("insertion") || otype.equals("placement") || otype.equals("refinement") || otype.equals("refinement_from_annotation"))
			{
				aa_flag = Boolean.parseBoolean(args[args.length - 5]);
				//initialize scoring matrix
	            _used_scoreTable = _nucleotide_nomenclature_scoreTable;
	            if(aa_flag) _used_scoreTable = _aminoacid_scoreTable;
			}
			
			if((otype.equals("refinement") || otype.equals("refinement_from_annotation") || otype.equals("insertion")) && informat.equals("fasta"))
			{
				if(!args[args.length - 4].equals("false")) checkpoint_prefix = args[args.length - 4];
			}
			
			outfn = args[args.length - 3]; //outputfile: a tree (insertion), jplace (placement), tsv (annotation/annotation_details)
			
			printDisInfoOnScreen = Boolean.parseBoolean(args[args.length - 2]); //is print to screen or not
			
			THREAD_NUM = Runtime.getRuntime().availableProcessors();
			int threads_N = Integer.parseInt(args[args.length - 1]);
			if(threads_N > 0) THREAD_NUM = threads_N;
			forkJoinPool = new ForkJoinPool(THREAD_NUM);
			String output_folder = getFolder(outfn);

			HashMap<Integer, String> queryList = null;
			HashMap<String, ArrayList<String>> clade2samples = null;

			// read tree
			NewickImporter tni = new NewickImporter(new FileReader(intfn));
			Tree tree = tni.importTree(null);
			TIPars myAdd = new TIPars(tree, otype, output_folder);
			Tree outtree = null;
						
			if (otype.equals("refinement") || otype.equals("refinement_from_annotation") || otype.equals("insertion") || otype.equals("placement"))
			{
				/////// read vcf file
				if (informat.contains("vcf") || informat.contains("Vcf") || informat.contains("VCF")) {
					multationSequencesMap.clear();
					seqIdxMap.clear();
					Arrays.fill(ref_sequence, (byte) '\0');
					readVCFAlignmentFile(insfn);
					readVCFAlignmentFile(inafn);
					if (otype.equals("placement") || otype.equals("insertion"))
						queryList = readVCFFile2Alignment(inqfn);
				} else {
					//// Read fasta file
					stringSequencesList.clear();
					seqIdxMap.clear();
					/// input msa
					readFastaAlignmentFile(insfn);
					/// input ancestral
					readFastaAlignmentFile(inafn);
					/// input query
					if (otype.equals("placement") || otype.equals("insertion"))
						queryList = readFastaFile2Alignment(inqfn);
				}
			}
			HashSet<String> taxaNames = null;		
			if (otype.equals("annotation") || otype.equals("annotation_details") || otype.equals("refinement") || otype.equals("refinement_from_annotation") || otype.equals("tree_BFS"))
			{
				taxaNames = setupHashSetOfTaxaName(tree);
				clade2samples = readCladeSample(inqfn, taxaNames);
			}
			
			HashMap<String, FlexibleNode> seqName2node = null;
			HashMap<FlexibleNode, String> assignmentNode2Clade = null;
			if (otype.equals("annotation_details") || otype.equals("tree_BFS"))
			{
				seqName2node = setupHashtableOfseqName2node(tree);
				if(otype.equals("annotation_details"))
					assignmentNode2Clade = readAnnotations(annotation_assignment, seqName2node);
			}
			
			HashMap<String, String> assignmentNodeName2Clade = null;
			if (otype.equals("refinement_from_annotation"))
			{
				assignmentNodeName2Clade = readAnnotations(annotation_assignment);
			}
			
			long startTime2 = System.currentTimeMillis();
			
			System.out.println("TIPars Version 2.0.0: Ultrafast tree annotation for Pango lineages of SRAS-CoV-2 sequences and application for tree refinement");
			System.out.println("Progress: " + otype);

			String progressInfo = "";
			// mutiple taxa insertion for vcf input
			if (otype.equals("insertion") && informat.contains("vcf")) {
				System.out.println("TreeFile: " + intfn);
	            System.out.println("TaxaFile: " + insfn);
	            System.out.println("AncestralSequenceFile: " + inafn);
	            System.out.println("QueryFile: " + inqfn);

				// init TIPars
				ArrayList<Integer> queryIdxs = new ArrayList<Integer>(queryList.keySet());
				queryIdxs.sort(Comparator.naturalOrder());
				for (int i = 0; i < queryIdxs.size(); i++) {
					int idx = queryIdxs.get(i);
					String name = queryList.get(idx);
					ConcurrentHashMap<Integer, Byte> query = multationSequencesMap.get(idx);
					String qid = "q" + (i + 1);
					String pid = "p" + (i + 1);
					outtree = myAdd.addQuerySequence(name, query, qid, pid, printDisInfoOnScreen, new double[3], otype, 0);
					myAdd.mytree = outtree;
					myAdd.setupHashtableOfnode2seqName();

					if (!printDisInfoOnScreen) {
						for (int j = 0; j < progressInfo.length(); j++)
							System.out.print("\b");
						progressInfo = (i + 1) + "/" + queryIdxs.size();
						System.out.print(progressInfo);
						if (i == queryIdxs.size() - 1)
							System.out.println();
					}
				}
			}
			// mutiple taxa insertion for fasta input
			else if (otype.equals("insertion")) {
				System.out.println("TreeFile: " + intfn);
	            System.out.println("TaxaFile: " + insfn);
	            System.out.println("AncestralSequenceFile: " + inafn);
	            System.out.println("QueryFile: " + inqfn);
	            
	            ArrayList<Integer> queryIdxs = new ArrayList<Integer>(queryList.keySet());
				queryIdxs.sort(Comparator.naturalOrder());
				
				int current_max_used_index4internalnode = -1;
				if(checkpoint_prefix != "")
				{
					String inFASTAFilePath = output_folder + "/" + checkpoint_prefix + "_addedInternodeP.fasta";
					File checkfile = new File(inFASTAFilePath);
					if (checkfile.exists()) {
						
						if(informat.contains("fasta"))
						{
							readFastaAlignmentFile(inFASTAFilePath);
						}
						
						System.out.println("AddedAncestralSequenceFile: " + inFASTAFilePath);
						
						String inTreePath = output_folder  + "/" + checkpoint_prefix + "_addedTree.nwk";
						tni = new NewickImporter(new FileReader(inTreePath));
						tree = tni.importTree(null);
						myAdd = new TIPars(tree, otype, output_folder);

						String inFilePath = output_folder + "/" + checkpoint_prefix + "_remainingSeqName.txt";
						ArrayList<String> innames = readSeqName(inFilePath);
						HashSet<String> innames_set = new HashSet<String>(innames);
						ArrayList<Integer> temp = new ArrayList<Integer>();
						for(int i = 0; i < queryIdxs.size(); i++)
						{
							int idx = queryIdxs.get(i);
							String name = queryList.get(idx);
							if(innames_set.contains(name)) temp.add(queryIdxs.get(i));
						}
						queryIdxs = temp;
					
						if(queryIdxs.size() != innames.size())
						{
							System.out.println("Error: reading checkpint query seq does not match the tree!");
						}

						for(HashMap.Entry<String, Integer> entry : seqIdxMap.entrySet()) {
							if(entry.getKey().charAt(0) == 'p')
							{
								int usedIdx = Integer.parseInt(entry.getKey().substring(1));
								if(usedIdx > current_max_used_index4internalnode) current_max_used_index4internalnode = usedIdx;
							}
						}
						current_max_used_index4internalnode = current_max_used_index4internalnode;
					}
					else
					{
						System.out.println("Warning: no checkpoint tree file and start from scratch!");
					}
				}
				
				for (int i = 0; i < queryIdxs.size(); i++) {
					int idx = queryIdxs.get(i);
					String name = queryList.get(idx);
					String query = stringSequencesList.get(idx);
					current_max_used_index4internalnode++;
					String qid = "q" + current_max_used_index4internalnode;
					String pid = "p" + current_max_used_index4internalnode;
					outtree = myAdd.addQuerySequence(name, query, qid, pid, printDisInfoOnScreen, new double[3], otype, 0);
					myAdd.mytree = outtree;
					myAdd.setupHashtableOfnode2seqName();

					if (!printDisInfoOnScreen) {
						for (int j = 0; j < progressInfo.length(); j++)
							System.out.print("\b");
						progressInfo = (current_max_used_index4internalnode + 1) + "/" + queryIdxs.size();
						System.out.print(progressInfo);
						if (i == queryIdxs.size() - 1)
							System.out.println();
					}
					
					//checkpoint to store the current outtree and remaining samples to be added
					if(checkpoint_prefix != "" && (i % 1000 == 999 || i == queryIdxs.size() - 1))
					{
						String outTreePath = output_folder  + "/" + checkpoint_prefix + "_addedTree.nwk";
						writeToTree(outtree, outTreePath, "insertion");
						String outFilePath = output_folder + "/" + checkpoint_prefix + "_remainingSeqName.txt";
						ArrayList<String> remainingSeqName = new ArrayList<String>();
						for(int j = i+1; j < queryIdxs.size(); ++j)
						{
							idx = queryIdxs.get(j);
							name = queryList.get(idx);
							remainingSeqName.add(name);
						}
						writeSeqName(outFilePath, remainingSeqName);
						
						if(informat.contains("fasta"))
						{
							String outFASTAFilePath = output_folder + "/" + checkpoint_prefix + "_addedInternodeP.fasta";
							ArrayList<String> desc = new ArrayList<String>();
							ArrayList<String> seq = new ArrayList<String>();
							for(HashMap.Entry<String, Integer> entry : seqIdxMap.entrySet()) {
								if(entry.getKey().charAt(0) == 'p')
								{
									desc.add(entry.getKey());
									seq.add(stringSequencesList.get(seqIdxMap.get(entry.getKey())));
								}
							}
							writeFASTA(desc, seq, outFASTAFilePath);
							System.out.println("addedInternodeP number: " + desc.size());
						}
					}
				}
			}
			// mutiple taxa placement for vcf input
			else if (otype.equals("placement") && informat.contains("vcf")) {
				System.out.println("TreeFile: " + intfn);
	            System.out.println("TaxaFile: " + insfn);
	            System.out.println("AncestralSequenceFile: " + inafn);
	            System.out.println("QueryFile: " + inqfn);
	            
				ArrayList<Integer> queryIdxs = new ArrayList<Integer>(queryList.keySet());
				queryIdxs.sort(Comparator.naturalOrder());
				placements = new String[queryIdxs.size()];
				for (int i = 0; i < queryIdxs.size(); i++) {
					int idx = queryIdxs.get(i);
					String name = queryList.get(idx);
					ConcurrentHashMap<Integer, Byte> query = multationSequencesMap.get(idx);
					String qid = "q" + (i + 1);
					String pid = "p" + (i + 1);
					outtree = myAdd.addQuerySequence(name, query, qid, pid, printDisInfoOnScreen, new double[3], otype, i);

					if (!printDisInfoOnScreen) {
						for (int j = 0; j < progressInfo.length(); j++)
							System.out.print("\b");
						progressInfo = (i + 1) + "/" + queryIdxs.size();
						System.out.print(progressInfo);
						if (i == queryIdxs.size() - 1)
							System.out.println();
					}
				}
			}
			// mutiple taxa placement for fasta input
			else if (otype.equals("placement")) {
				System.out.println("TreeFile: " + intfn);
	            System.out.println("TaxaFile: " + insfn);
	            System.out.println("AncestralSequenceFile: " + inafn);
	            System.out.println("QueryFile: " + inqfn);
	            
				ArrayList<Integer> queryIdxs = new ArrayList<Integer>(queryList.keySet());
				queryIdxs.sort(Comparator.naturalOrder());
				placements = new String[queryIdxs.size()];
				for (int i = 0; i < queryIdxs.size(); i++) {
					int idx = queryIdxs.get(i);
					String name = queryList.get(idx);
					String query = stringSequencesList.get(idx);
					String qid = "q" + (i + 1);
					String pid = "p" + (i + 1);
					outtree = myAdd.addQuerySequence(name, query, qid, pid, printDisInfoOnScreen, new double[3], otype, i);

					if (!printDisInfoOnScreen) {
						for (int j = 0; j < progressInfo.length(); j++)
							System.out.print("\b");
						progressInfo = (i + 1) + "/" + queryIdxs.size();
						System.out.print(progressInfo);
						if (i == queryIdxs.size() - 1)
							System.out.println();
					}
				}
			} 
			else if (otype.equals("annotation")) 
			{
				System.out.println("TreeFile: " + intfn);
	            System.out.println("LabelFile: " + inqfn);
	            System.out.println("F1score threshold: " + fscore_min);
	            
				HashMap<String, String> assignmentMap = myAdd.treeAnnotation(clade2samples, fscore_min, printDisInfoOnScreen);
				writeAssignClade(outfn, assignmentMap);
				System.out.println("INFO: write annotation assignment to file " + outfn);
			} 
			else if (otype.equals("annotation_details")) 
			{
				System.out.println("TreeFile: " + intfn);
	            System.out.println("LabelFile: " + inqfn);
	            System.out.println("AnnotationFile: " + annotation_assignment);

				HashMap<FlexibleNode, ScoreAndStringList> assignedNode2Samples = myAdd.annotationStatistics(clade2samples, assignmentNode2Clade, seqName2node, taxaNames, printDisInfoOnScreen, true);
				writeAssignCladeAndTaxa(outfn, assignmentNode2Clade, assignedNode2Samples);
				System.out.println("INFO: write annotation details to file " + outfn);
				
				final String JSON_OUTPUT_FILE_PATH = output_folder.isEmpty() ? "graph-data-generated.js" : output_folder + "/graph-data-generated.js"; //changeFileExtension(outfn, "js");
				final String ANNOTATION_TOOLTIP_NAME = "Pango lineage";
				generateGraphJSON(
					JSON_OUTPUT_FILE_PATH, 
					ANNOTATION_TOOLTIP_NAME,
					outfn,
					clade2samples
				);
				copyStaticWebFiles(output_folder);
				System.out.println("INFO: write annotation visualization (without unannotated samples) to file " + JSON_OUTPUT_FILE_PATH);
				System.out.println("INFO: please click the 'graph.html' in the output folder for visualization." );
			}
			else if (otype.equals("refinement"))
			{
				System.out.println("TreeFile: " + intfn);
				System.out.println("TaxaFile: " + insfn);
	            System.out.println("LabelFile: " + inqfn);
	            
				outtree = myAdd.doTreeRefinement(clade2samples, informat, printDisInfoOnScreen, null, checkpoint_prefix);
			}
			else if (otype.equals("refinement_from_annotation"))
			{
				System.out.println("TreeFile: " + intfn);
				System.out.println("TaxaFile: " + insfn);
	            System.out.println("LabelFile: " + inqfn);
	            System.out.println("AnnotationFile: " + annotation_assignment);
				outtree = myAdd.doTreeRefinement(clade2samples, informat, printDisInfoOnScreen, assignmentNodeName2Clade, checkpoint_prefix);
			}
			else if (otype.equals("tree_BFS"))
			{
				System.out.println("TreeFile: " + intfn);
	            System.out.println("LabelFile: " + inqfn);
	            System.out.println("ExploreTreeNodeLimit: " + exploreTreeNodeLimit);
	            System.out.println("SmallBubbleLimit: " + smallBubbleLimit);
	            System.out.println("SmallClusterLimit: " + smallClusterLimit);
	            System.setProperty("java.util.Arrays.useLegacyMergeSort", "true");
	            HashMap<String, String> assignmentMap = myAdd.treeAnnotation(clade2samples, fscore_min, printDisInfoOnScreen);
	            assignmentNode2Clade = new HashMap<FlexibleNode, String>();
	            for (Entry<String,String> entry : assignmentMap.entrySet()) {
	            	if(seqName2node.containsKey(entry.getKey()))
	            		assignmentNode2Clade.put(seqName2node.get(entry.getKey()), entry.getValue());
	            	else
	            		System.out.println("WARNNING: can not find assigned node " + entry.getKey());
	            }
	            
	            HashMap<FlexibleNode, ScoreAndStringList> assignedNode2Samples = myAdd.annotationStatistics(clade2samples, assignmentNode2Clade, seqName2node, taxaNames, printDisInfoOnScreen, false);
	            HashMap<FlexibleNode, HashMap<String, HashSet<String>>> stratification = stratifyTree(tree, assignedNode2Samples, exploreTreeNodeLimit, smallBubbleLimit, smallClusterLimit, seqName2node);
				writeAllStratification2TSV(stratification, assignedNode2Samples, assignmentNode2Clade, seqName2node, outfn);
			}
			
			long endTime2 = System.currentTimeMillis();
			long totalTime2 = endTime2 - startTime2;
			
			//write nwk tree or jplace file
			if (outtree != null) 
			{
				writeToTree(outtree, outfn, otype);
				System.out.println("INFO: output tree/placement to file " + outfn);
			}

			long endTime = System.currentTimeMillis();
			long totalTime = endTime - startTime;

			System.out.println("Processing time (excluding I/O, seconds): " + (double) totalTime2 / 1000);
			System.out.println("Overall time (seconds): " + (double) totalTime / 1000);
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

