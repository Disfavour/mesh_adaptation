// Tree.h -- left red-black tree
//

#pragma once


#ifdef _DEBUG

inline void assert(bool test)
{
	if(!test)__asm int 3;
}

#else

inline void assert(bool test)
{
}

#endif



class NoCopy
{
	NoCopy(NoCopy &);
	NoCopy &operator = (NoCopy &);

public:
	NoCopy()
	{
	}
};


class BaseTreeNode : public NoCopy
{
	enum NodeType
	{
		t_root = 0, t_right, t_left, t_left_red
	};


	BaseTreeNode *parent, *left, *right;
	NodeType type;


	BaseTreeNode *First()
	{
		for(BaseTreeNode *node = this;; node = node->left)if(!node->left)return node;
	}

	BaseTreeNode *Last()
	{
		for(BaseTreeNode *node = this;; node = node->right)if(!node->right)return node;
	}


	void SetParent(BaseTreeNode *ptr, NodeType tp)
	{
		parent = ptr;  type = tp;
		switch(tp)
		{
		case t_root:  *reinterpret_cast<BaseTreeNode **>(ptr) = this;  return;
		case t_right:  ptr->right = this;  return;
		default:  ptr->left = this;  return;
		}
	}

	void RotateLeft()
	{
		BaseTreeNode *next = right;  next->SetParent(parent, type);
		if(right = next->left)
		{
			right->parent = this;  right->type = t_right;
		}
		next->left = this;  parent = next;  type = t_left_red;
	}

	void RotateRight()
	{
		BaseTreeNode *next = left;  next->SetParent(parent, type);
		if(left = next->right)
		{
			left->parent = this;  left->type = t_left;
		}
		next->right = this;  parent = next;  type = t_right;
	}

	BaseTreeNode *MakeRed()
	{
		if(type == t_left)
			if(parent->type == t_left_red)
			{
				parent->parent->RotateRight();  return parent;
			}
			else type = t_left_red;
		else if(type == t_right)
			if(parent->left && parent->left->type == t_left_red)
			{
				parent->left->type = t_left;  return parent;
			}
			else
			{
				parent->RotateLeft();
				if(type == t_left_red)
				{
					left->type = t_left;  parent->RotateRight();  return this;
				}
			}
		return 0;
	}

	BaseTreeNode *RemoveRed()
	{
		if(type == t_left_red)type = t_left;
		else if(type == t_left)
			if(parent->right->left && parent->right->left->type == t_left_red)
			{
				parent->right->RotateRight();  parent->RotateLeft();  parent->type = t_left;
			}
			else
			{
				parent->RotateLeft();  return parent->parent;
			}
		else if(type == t_right)
			if(parent->left->type == t_left)
				if(parent->left->left && parent->left->left->type == t_left_red)
				{
					parent->left->left->type = t_left;  parent->RotateRight();
				}
				else
				{
					parent->left->type = t_left_red;  return parent;
				}
			else
				if(parent->left->right->left && parent->left->right->left->type == t_left_red)
				{
					parent->left->RotateLeft();  parent->RotateRight();
				}
				else
				{
					parent->RotateRight();  parent->left->type = t_left_red;
				}
		return 0;
	}

	void Swap(BaseTreeNode *node)
	{
		BaseTreeNode *l = left, *r = right, *p = parent;  NodeType t = type;
		if(node->parent == this)
		{
			parent = node;  type = node->type;  node->SetParent(p, t);
			if(type == t_right)
			{
				if(left = node->left)left->parent = this;  if(node->left = l)l->parent = node;
				if(right = node->right)right->parent = this;  node->right = this;
			}
			else
			{
				if(left = node->left)left->parent = this;  node->left = this;
				if(right = node->right)right->parent = this;  if(node->right = r)r->parent = node;
			}
		}
		else
		{
			SetParent(node->parent, node->type);  node->SetParent(p, t);
			if(left = node->left)left->parent = this;  if(node->left = l)l->parent = node;
			if(right = node->right)right->parent = this;  if(node->right = r)r->parent = node;
		}
	}

	/*void Check()  // DEBUG
	{
		if(type == t_left_red && parent->type == t_left_red)__asm int 3;
		if(left)
		{
			if(left->parent != this || left->type <= t_right)__asm int 3;
			if(!right && (left->type == t_left || left->left || left->right))__asm int 3;
			left->Check();
		}
		if(right)
		{
			if(right->parent != this || right->type != t_right || !left)__asm int 3;
			right->Check();
		}
		if(type == t_root && *reinterpret_cast<BaseTreeNode **>(parent) != this)__asm int 3;
	}*/


public:
	template<class T> struct Deleter
	{
		void operator () (T *del)
		{
			del->parent = 0;  delete del;
		}
	};

	template<class T> struct Remover
	{
		void operator () (T *del)
		{
			del->parent = 0;
		}
	};

	class Pos
	{
		BaseTreeNode *node;
		bool after;

		Pos() : node(0)
		{
		}

	public:
		Pos(BaseTreeNode *ptr, bool aft)
		{
			if(!ptr)node = 0;
			else if(aft)
				if(ptr->right)
				{
					node = ptr->right->First();  after = false;
				}
				else
				{
					node = ptr;  after = true;
				}
			else
				if(ptr->left)
				{
					node = ptr->left->Last();  after = true;
				}
				else
				{
					node = ptr;  after = false;
				}
		}

		Pos(Pos &pos) : node(pos.node), after(pos.after)
		{
		}

		Pos &operator = (Pos &pos)
		{
			node = pos.node;  after = pos.after;  return *this;
		}

		BaseTreeNode *Node()
		{
			return node;
		}

		bool After()
		{
			return after;
		}

		friend class BaseTree;
	};


	BaseTreeNode() : parent(0)
	{
	}

	~BaseTreeNode()
	{
		if(parent)Remove();
	}

	BaseTreeNode *Prev()
	{
		assert(parent != 0);  BaseTreeNode *node = this;
		if(!left)
		{
			while(node->type > t_right)node = node->parent;
			return node->type ? node->parent : 0;
		}
		for(node = node->left;; node = node->right)if(!node->right)return node;
	}

	BaseTreeNode *Next()
	{
		assert(parent != 0);  BaseTreeNode *node = this;
		if(!right)
		{
			while(node->type == t_right)node = node->parent;
			return node->type ? node->parent : 0;
		}
		for(node = node->right;; node = node->left)if(!node->left)return node;
	}


	void Insert(BaseTreeNode *node, bool after)
	{
		assert(node->parent == 0);  node->parent = this;  node->left = node->right = 0;
		if(after)
		{
			assert(right == 0);  right = node;  node->type = t_right;
		}
		else
		{
			assert(left == 0);  left = node;  node->type = t_left;
		}
		do node = node->MakeRed();
		while(node);
	}

	void Remove()
	{
		assert(parent != 0);  while(left)Swap(left->Last());
		
		BaseTreeNode *node = this;
		do node = node->RemoveRed();
		while(node);

		if(type > t_right)parent->left = 0;
		else if(type == t_right)parent->right = 0;
		else *reinterpret_cast<BaseTreeNode **>(parent) = 0;
		parent = 0;
	}

	bool Used()
	{
		return parent != 0;
	}


	friend class BaseTree;
};

class BaseTree : public NoCopy
{
	BaseTreeNode *root;


public:
	typedef BaseTreeNode::Pos Pos;


	/*void Check()  // DEBUG
	{
		if(root)
		{
			root->Check();  if(root->type != BaseTreeNode::t_root)__asm int 3;
		}
	}*/


	BaseTree() : root(0)
	{
	}

	BaseTreeNode *First()
	{
		return root ? root->First() : 0;
	}

	BaseTreeNode *Last()
	{
		return root ? root->Last() : 0;
	}


	template<class T, typename D> void Clear(D del)
	{
		if(root)for(BaseTreeNode *node = root;;)
		{
			if(node->left)node = node->left;
			else if(node->right)node = node->right;
			else
			{
				BaseTreeNode *ptr = node;  node = node->parent;
				if(ptr->type > BaseTreeNode::t_right)node->left = 0;
				else if(ptr->type == BaseTreeNode::t_right)node->right = 0;
				else
				{
					del(static_cast<T *>(ptr));  return;
				}
				del(static_cast<T *>(ptr));
			}
		}
	}

	template<class T, typename C> Pos Find(C cmp)
	{
		Pos res;
		for(BaseTreeNode *node = root; node;)
			if(cmp(static_cast<T *>(node)))
			{
				res.node = node;  res.after = true;  node = node->right;
			}
			else
			{
				res.node = node;  res.after = false;  node = node->left;
			}
		return res;
	}

	void Insert(BaseTreeNode *node, Pos pos)
	{
		if(pos.node)pos.node->Insert(node, pos.after);
		else
		{
			root = node;  node->parent = reinterpret_cast<BaseTreeNode *>(&root);
			node->left = node->right = 0;  node->type = BaseTreeNode::t_root;
		}
	}

	bool operator ! ()
	{
		return root == 0;
	}

	operator bool()
	{
		return root != 0;
	}
};


template<class T> class TreeNode : public BaseTreeNode
{
public:
	T *Prev()
	{
		return static_cast<T *>(BaseTreeNode::Prev());
	}

	T *Next()
	{
		return static_cast<T *>(BaseTreeNode::Next());
	}
};

template<class T, typename D = BaseTreeNode::Remover<T> > class Tree : public BaseTree
{
	D del;


public:
	class Pos : public BaseTree::Pos
	{
	public:
		Pos(T *node, bool after) : BaseTree::Pos(node, after)
		{
		}

		Pos(BaseTree::Pos &pos) : BaseTree::Pos(pos)
		{
		}

		T *Node()
		{
			return static_cast<T *>(BaseTree::Pos::Node());
		}
	};

	struct BaseCmp
	{
		T *node;

		BaseCmp(T *ptr) : node(ptr)
		{
		}

		bool operator () (T *cmp)
		{
			return *cmp < *node;
		}
	};


	Tree()
	{
	}

	Tree(D &d) : del(d)
	{
	}

	~Tree()
	{
		Clear<T>(del);
	}

	T *First()
	{
		return static_cast<T *>(BaseTree::First());
	}

	T *Last()
	{
		return static_cast<T *>(BaseTree::Last());
	}

	template<typename C> Pos Find(C cmp)
	{
		return BaseTree::Find<T, C>(cmp);
	}

	void Insert(T *node, Pos pos)
	{
		BaseTree::Insert(node, pos);
	}

	void Insert(T *node)
	{
		BaseTree::Insert(node, BaseTree::Find<T>(BaseCmp(node)));
	}
};

template<class T> class OwningTree : public Tree<T, BaseTreeNode::Deleter<T>>
{
};