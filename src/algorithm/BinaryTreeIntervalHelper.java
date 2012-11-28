package algorithm;

import dataStructure.IntervalAbs;

public class BinaryTreeIntervalHelper {
	// Root node pointer. Will be null for an empty tree.
	private Node root;

	/*
	 * --Node-- The binary tree is built using this nested node class. Each node
	 * stores one data element, and has left and right sub-tree pointer which
	 * may be null. The node is a "dumb" nested class -- we just use it for
	 * storage; it does not have any methods.
	 */
	private static class Node {
		Node left;
		Node right;
		IntervalAbs data;

		Node(IntervalAbs newData) {
			left = null;
			right = null;
			data = newData;
		}
	}

	/**
	 * Creates an empty binary tree -- a null root pointer.
	 */
	public BinaryTreeIntervalHelper() {
		root = null;
	}

	/**
	 * Returns true if the given target is in the binary tree. Uses a recursive
	 * helper.
	 */
	public boolean lookup(IntervalAbs data) {
		return (lookup(root, data));
	}

	/**
	 * Recursive lookup -- given a node, recur down searching for the given
	 * data.
	 */
	private boolean lookup(Node node, IntervalAbs data) {
		if (node == null) {
			return (false);
		}

		if (data == node.data) {
			return (true);
		} else if (data.getStartAbs() < node.data.getStartAbs()) {
			return (lookup(node.left, data));
		} else {
			return (lookup(node.right, data));
		}
	}

	/**
	 * Inserts the given data into the binary tree. Uses a recursive helper.
	 */
	public void insert(IntervalAbs data) {
		root = insert(root, data);
	}

	/**
	 * Recursive insert -- given a node pointer, recur down and insert the given
	 * data into the tree. Returns the new node pointer (the standard way to
	 * communicate a changed pointer back to the caller).
	 */
	private Node insert(Node node, IntervalAbs data) {
		if (node == null) {
			node = new Node(data);
		} else {
			if (data.getStartAbs() <= node.data.getStartAbs()) {
				node.left = insert(node.left, data);
			} else {
				node.right = insert(node.right, data);
			}
		}

		return (node); // in any case, return the new pointer to the caller
	}

	/**
	 * Returns the number of nodes in the tree. Uses a recursive helper that
	 * recurs down the tree and counts the nodes.
	 */

	public int size() {
		return (size(root));
	}

	private int size(Node node) {
		if (node == null)
			return 0;
		else
			return (size(node.left) + 1 + size(node.right));

	}

	/**
	 * Returns the max root-to-leaf depth of the tree. Uses a recursive helper
	 * that recurs down to find the max depth.
	 */
	public int maxDepth() {
		return (maxDepth(root));
	}

	private int maxDepth(Node node) {
		if (node == null)
			return 0;
		else {
			int lDepth = maxDepth(node.left);
			int rDepth = maxDepth(node.right);

			// use the larger + 1
			return (Math.max(lDepth, rDepth) + 1);
		}
	}

	/**
	 * Returns the min value in a non-empty binary search tree. Uses a helper
	 * method that iterates to the left to find the min value.
	 */
	public IntervalAbs minValue() {
		return (minValue(root));
	}

	/**
	 * Finds the min value in a non-empty binary search tree. the leftmost leaf
	 * node.data which node.left is null has the minValue.
	 */
	private IntervalAbs minValue(Node node) {
		Node current = node;
		while (current.left != null)
			current = current.left;
		return current.data;
	}

	/**
	 * Returns the max value in a non-empty binary search tree. Uses a helper
	 * method that iterates to the right to find the max value.
	 */
	public IntervalAbs maxValue() {
		return (maxValue(root));
	}

	/**
	 * Finds the max value in a non-empty binary search tree. the rightmost leaf
	 * node.data which node.right is null has the maxValue.
	 */
	private IntervalAbs maxValue(Node node) {
		Node current = node;
		while (current.right != null)
			current = current.right;
		return current.data;
	}

	/**
	 * Prints the node values in the "inorder" order. Uses a recursive helper to
	 * do the traversal.
	 */
	public void printTree() {
		printTree(root);
		System.out.println();
	}

	private void printTree(Node node) {
		if (node == null)
			return;

		// left, node itself, right
		printTree(node.left);
		System.out.print(node.data + "  ");
		printTree(node.right);
	}

	/**
	 * Prints the node values in the "postorder" order. Uses a recursive helper
	 * to do the traversal.
	 */
	public void printPostorder() {
		printPostorder(root);
		System.out.println();
	}

	public void printPostorder(Node node) {
		if (node == null)
			return;

		// first recur on both subtrees
		printPostorder(node.left);
		printPostorder(node.right);

		// then deal with the node
		System.out.print(node.data + "  ");
	}



	/**
	 * Given a binary tree, prints out all of its root-to-leaf paths, one per
	 * line. Uses a recursive helper to do the work.
	 */
	public void printPaths() {
		long[] path = new long[1000000];
		printPaths(root, path, 0);
		// for (int i : path)
		// System.out.println(i);
	}

	/**
	 * Recursive printPaths helper -- given a node, and an array containing the
	 * path from the root node up to but not including this node, prints out all
	 * the root-leaf paths.
	 */
	private void printPaths(Node node, long[] path, int pathLen) {
		if (node == null)
			return;
		path[pathLen] = node.data.getStartAbs();
		pathLen++;
		// it's a leaf, so print the path that led to here. It is a leaf, if
		// node.left AND node.roght = NULL
		if (node.left == null && node.right == null) {
			for (int i = 0; i < pathLen; i++) {
				System.out.print(path[i] + " ");
			}
			System.out.println();
		} else {
			// otherwise try both subtrees
			printPaths(node.left, path, pathLen);
			printPaths(node.right, path, pathLen);
		}
	}

	/**
	 * Changes the tree into its mirror image.
	 * 
	 * So the tree... 4 / \ 2 5 / \ 1 3
	 * 
	 * is changed to... 4 / \ 5 2 / \ 3 1
	 * 
	 * Uses a recursive helper that recurs over the tree, swapping the
	 * left/right pointers.
	 */
	public void mirror() {
		mirror(root);
	}

	private void mirror(Node node) {
		if (node != null) {
			// do the sub-trees
			mirror(node.left);
			mirror(node.right);
			Node temp = node.left;
			node.left = node.right;
			node.right = temp;
		}
	}

	public void printInOrder(Node node) {
		if (node != null) {
			printInOrder(node.left);
			System.out.println("  Traversed " + node.data);
			printInOrder(node.right);
		}
	}
}
