#pragma once
#include <vector>
#include <stack>

// author: Kaan Ucar

// Node of KDTree that represents a axis-aligned bounding box.
//
// Nodes eighter have all the children or do not have any children. But depth of leaf nodes might show variety.
// Since dimension is fixed once the node created min-max bounds and center coordinates stored as Scalar array.
// Dimension must be between in [1, bitsize of int - 1] interval. For example when sizeof int is 4 byte maximum dimension size supported is 31.
// Behaviour is undefined when dimension is out of interval.
//
// typename Scalar : floating point type for min-max bounds
// int Dim         : dimension[1, bitsize of int]
// typename Object : object type store
template<typename Scalar, int Dim, typename Object>
class KDTreeNode
{
public:
	// Constructor
	//
	// Warning! : Do not create nodes by constructor unless you creating a root. Use addChildren functions instead.
	//
	// min-max bounds must have Dim elements in them.
	// Pass parent as nullptr when it is a root node.
	// Once the node created min-max bounds and center of it cannot be changed.
	//
	// minBound : coordinates of min bound
	// maxBound : coordiantes of max bound
	// parent   : parent node
	// object   : object to point (optional)
	KDTreeNode(const Scalar *minBound, const Scalar *maxBound, KDTreeNode *parent, const Object& object);

	// Destructor
	// Removes all the children nodes
	virtual ~KDTreeNode();

	// return : dimension
	int dimension() const { return Dim; }

	// return : min bound coordinates
	const Scalar* getMinBounds() const { return m_minBound; }

	// return : max bound coordinates
	const Scalar* getMaxBounds() const { return m_maxBound; }

	// return : center coordinates
	const Scalar* getCenter() const { return m_center; }

	// return : true if the node has children.
	bool hasChildren() const { return m_children != nullptr; }

	// User responsible to check whether the node has children or not.
	//
	// return : child node at given index
	KDTreeNode* getChild(const int childIndex) const { return m_children[childIndex]; }

	// return : true if the node is not a root node
	bool hasParent() const { return m_parent != nullptr; }

	// return : parent
	KDTreeNode* getParent() const { return m_parent; }

	// Even if node does not have any children returns potential number of childrens.
	//
	// return : number of childrens(2^dimension)
	int numChildren() const { return 1 << Dim; }

	// Add children to node.
	//
	// Children's bound and center values calculated automatically.
	// childrenObjects must be ordered respect to Morton(Z) order or it can be passed as nullptr to initialize all the children object by default constructor.
	//
	// childrenObjects : list of object pointers.
	void addChildren(const Object *childrenObjects = nullptr);

	// Delete all the children nodes recursively.
	void removeChildren(const bool deleteObjects = false);

	// User responsible to check whether the node has parent or not.
	//
	// return : index in parent.
	int getIndex() const;

	// Get list of leaf nodes at given direction. Return itself if does not have any children.
	// For example getting leaf nodes that closest to -y direction, can be achieved by passing parameters as axis = 1, sign = false.
	//      x y z ...
	// axis 0 1 2 ...
	// 
	//      -     +
	// sign false true
	//
	// axis : pole axis index
	// sign : sign of axis, true when positive
	// return : list of pole leaf nodes
	std::vector<KDTreeNode*> getPoles(const int axis, const bool sign) const;

	// Get list of neighbour leaf nodes at given direction.
	// For example getting neigbour leaf nodes at +z direction, can be achieved by passing parameters as axis = 2, sign = true.
	//      x y z ...
	// axis 0 1 2 ...
	// 
	//      -     +
	// sign false true
	//
	// axis : neighbour axis index
	// sign : sign of axis, true when positive
	// return : list of neighbour leaf nodes
	std::vector<KDTreeNode*> getNeighbours(const int axis, const bool sign) const;

	// Get list of all neighbour leaf nodes.
	//
	// return : list of all neighbour leaf nodes
	std::vector<KDTreeNode*> getAllNeighbours() const;

	Object object;

private:
	Scalar m_minBound[Dim];
	Scalar m_maxBound[Dim];
	Scalar m_center[Dim];

	// nullptr when node is a root node
	KDTreeNode *m_parent;

	// nullptr when does not have children for space saving
	KDTreeNode **m_children;
};

// Constructor
//
// Warning! : Do not create nodes by constructor unless you really have to do it. Use KDTree constructor or addChildren functions instead.
//
// minBound : coordinates of min bound
// maxBound : coordiantes of max bound
// parent   : parent node
// object   : object to store
template<typename Scalar, int Dim, typename Object>
KDTreeNode<Scalar, Dim, Object>::KDTreeNode(const Scalar *minBound, const Scalar *maxBound, KDTreeNode *parent, const Object& object)
	: object(object), m_parent(parent), m_children(nullptr)
{
	std::memcpy(m_minBound, minBound, sizeof(Scalar) * Dim);
	std::memcpy(m_maxBound, maxBound, sizeof(Scalar) * Dim);

	for (int i = 0; i < Dim; ++i)
		m_center[i] = (minBound[i] + maxBound[i]) / 2;
}

// Destructor
//
// Removes children. 
template<typename Scalar, int Dim, typename Object>
inline KDTreeNode<Scalar, Dim, Object>::~KDTreeNode()
{
	if (hasChildren())
		removeChildren();
}

// Add children to node.
//
// childrenObjects : list of object pointers.
template<typename Scalar, int Dim, typename Object>
void KDTreeNode<Scalar, Dim, Object>::addChildren(const Object *childrenObjects)
{
	const int NUM_CHILDREN = numChildren();
	m_children = new KDTreeNode<Scalar, Dim, Object>*[NUM_CHILDREN];

	Scalar minBound[Dim];
	Scalar maxBound[Dim];

	// place children respect to Morton(Z) order.
	for (int i = 0; i < NUM_CHILDREN; ++i)
	{
		for (int j = 0; j < Dim; ++j)
		{
			if (i & (1 << j))
			{
				minBound[j] = m_center[j];
				maxBound[j] = m_maxBound[j];
			}
			else
			{
				minBound[j] = m_minBound[j];
				maxBound[j] = m_center[j];
			}
		}

		if (childrenObjects != nullptr)
			m_children[i] = new KDTreeNode<Scalar, Dim, Object>(minBound, maxBound, this, childrenObjects[i]);
		else
			m_children[i] = new KDTreeNode<Scalar, Dim, Object>(minBound, maxBound, this, Object());
	}
}

// Remove all the children nodes recursively.
template<typename Scalar, int Dim, typename Object>
void KDTreeNode<Scalar, Dim, Object>::removeChildren(const bool deleteObjects)
{
	const int NUM_CHILDREN = numChildren();

	for (int i = 0; i < NUM_CHILDREN; ++i)
		delete m_children[i];

	delete m_children;
	m_children = nullptr;
}

// return : index in parent.
template<typename Scalar, int Dim, typename Object>
int KDTreeNode<Scalar, Dim, Object>::getIndex() const
{
	int index = 0;

	// set pozitif axis
	for (int i = 0; i < Dim; ++i)
		if (m_minBound[i] == m_parent->m_center[i])
			index |= (1 << i);

	return index;
}

// Get list of leaf nodes at given direction. Return itself if does not have any children.
//
// axis : pole axis index
// sign : sign of axis, true when positive
// return : list of pole leaf nodes
template<typename Scalar, int Dim, typename Object>
std::vector < KDTreeNode< Scalar, Dim, Object >* > KDTreeNode<Scalar, Dim, Object>::getPoles(const int axis, const bool sign) const
{
	std::vector<KDTreeNode<Scalar, Dim, Object>*> poles;

	// if no children return itself
	if (!hasChildren())
	{
		poles.push_back(const_cast<KDTreeNode<Scalar, Dim, Object>*>(this));
		return poles;
	}

	const int NUM_CHILDREN = numChildren();
	const int nthBit = 1 << axis;
	const int nthBitSign = sign ? nthBit : 0;

	for (int i = 0; i < NUM_CHILDREN; ++i)
		if (nthBitSign == (i & nthBit))
		{
			const auto& subPoles = m_children[i]->getPoles(axis, sign);
			poles.insert(std::end(poles), std::begin(subPoles), std::end(subPoles));
		}

	return poles;
}

// Get list of neighbour leaf nodes at given direction.
//
// axis : neighbour axis index
// sign : sign of axis, true when positive
// return : list of neighbour leaf nodes
template<typename Scalar, int Dim, typename Object>
std::vector < KDTreeNode< Scalar, Dim, Object >* > KDTreeNode<Scalar, Dim, Object>::getNeighbours(const int axis, const bool sign) const
{
	const int nthBit = 1 << axis;
	const int nthBitSign = sign ? nthBit : 0;

	std::stack<int> indexes;
	auto node = this;

	do {
		if (!node->hasParent()) // no neighbour found, node should be at the edge of kdtree
			return std::vector<KDTreeNode<Scalar, Dim, Object>*>();

		indexes.push(node->getIndex());
		node = node->m_parent;
	} while (nthBitSign == (indexes.top() & nthBit)); // loop until find the common parent

	while (!indexes.empty() && node->hasChildren())
	{
		node = node->getChild(indexes.top() ^ nthBit); // get neighbour node that has same or bigger size
		indexes.pop();
	}

	return node->getPoles(axis, !sign); // get closest leaf nodes
}

// Get list of all neighbour leaf nodes.
//
// return : list of all neighbour leaf nodes
template<typename Scalar, int Dim, typename Object>
inline std::vector<KDTreeNode<Scalar, Dim, Object>*> KDTreeNode<Scalar, Dim, Object>::getAllNeighbours() const
{
	std::vector<KDTreeNode<Scalar, Dim, Object>*> neighbours;

	// append neighbours at all direction
	for (int i = 0; i < Dim; ++i)
	{
		const auto& negativeNeighbours = getNeighbours(i, false);
		neighbours.insert(std::end(negativeNeighbours), std::begin(negativeNeighbours), std::end(negativeNeighbours));
		const auto& positiveNeighbours = getNeighbours(i, true);
		neighbours.insert(std::end(positiveNeighbours), std::begin(positiveNeighbours), std::end(positiveNeighbours));
	}

	return neighbours;
}
