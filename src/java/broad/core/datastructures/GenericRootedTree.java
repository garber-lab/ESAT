package broad.core.datastructures;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class GenericRootedTree<T> {
	private String name;
	private TreeNode<T> root;
	private ArrayList<ArrayList<TreeNode<T>>> levels;
	
	public GenericRootedTree(String name, T rootElement) {
		levels = new ArrayList<ArrayList<TreeNode<T>>>();
		try {
			this.root = new TreeNode<T>(this, rootElement);
		} catch (IllegalTreeOperationException e) {
			e.printStackTrace();
		}
		this.name = name;		
	}
	
	public TreeNode<T> getRoot() { return root;}
	public String getName() { return name;}
	
	private void addNodeToLevel(int level, TreeNode<T> node) throws IllegalTreeOperationException {
		if(level > levels.size()) {
			throw new IllegalTreeOperationException("Trying to add a level that is more than one level deeper than any existing  level ");
		}
		if(level == levels.size()) {
			ArrayList<TreeNode<T>> levelList = new ArrayList<TreeNode<T>>();
			levels.add(levelList);
		}
		ArrayList<TreeNode<T>> levelNodeList = levels.get(level);
		levelNodeList.add(level, node);
	}
	
	public List<TreeNode<T>> getNodes(int level) { return levels.get(level);}
	
	public List<TreeNode<T>> getLeaves() { return getNodes(levels.size() - 1); }
	
	public List<T> getNodeElements(int level) {
		List<TreeNode<T>> nodes = getNodes(level);
		ArrayList<T> elements = new ArrayList<T>(nodes.size());
		
		Iterator<TreeNode<T>> it = nodes.iterator();
		while(it.hasNext()) {
			elements.add(it.next().getElement());
		}
		return elements;
	}
	
	public static class TreeNode<T> {
		private TreeNode<T> parent;
		private ArrayList<TreeNode<T>> children;
		private T nodeElement;
		private int level;
		
		public TreeNode(GenericRootedTree<T> tree, T element) throws IllegalTreeOperationException {
			nodeElement = element;
			this.level = 0;
			children = new ArrayList<TreeNode<T>>();
			tree.addNodeToLevel(level, this);
		}
		public TreeNode(GenericRootedTree<T> tree, T element, TreeNode<T> parent) throws IllegalTreeOperationException {
			nodeElement = element;
			this.parent = parent;
			this.level = parent.getLevel() + 1;
			children = new ArrayList<TreeNode<T>>();
			tree.addNodeToLevel(level, this);
		}
		
		public int getLevel() { return level;}
		
		protected void setLevel(int level) { this.level = level;}
		
		public void setParent(TreeNode<T> parent) {
			this.parent = parent;
			parent.addChild(this);
			this.level = parent.getLevel() + 1;
			Iterator<TreeNode<T>> childIt = children.iterator();
			while(childIt.hasNext()) {
				TreeNode<T> child = childIt.next();
				child.setParent(this);
			}
		}

		public void addChild(TreeNode<T> child) {
			children.add(child);
		}
		
		public T getElement() { return nodeElement;}
		
		public TreeNode<T> getParent() { return this.parent;}
	}
}
