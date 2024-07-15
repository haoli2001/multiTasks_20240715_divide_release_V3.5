#include "tree2vector.h"
#include <string.h>
#include <stdlib.h>
#include <memory>



void _give_tree_index(KD_Node *node, int *index)
{
	if (node == NULL)
		return;
	node->index = *index;
	(*index)++;
	_give_tree_index(node->LeftChild, index);
	_give_tree_index(node->RightChild, index);
}

void _bianli_sz2(KD_Node *root, KD_Node_V *nodeg)
{
	if (root == NULL)
	{
		return;
	}
	else
	{
		root->LeftChild == NULL ? nodeg->LeftIndex = -1 : nodeg->LeftIndex = root->LeftChild->index;
		root->RightChild == NULL ? nodeg->RightIndex = -1 : nodeg->RightIndex = root->RightChild->index;
		root->rope[0] == NULL ? nodeg->RopeIndex[0] = -1 : nodeg->RopeIndex[0] = root->rope[0]->index;
		root->rope[1] == NULL ? nodeg->RopeIndex[1] = -1 : nodeg->RopeIndex[1] = root->rope[1]->index;
		root->rope[2] == NULL ? nodeg->RopeIndex[2] = -1 : nodeg->RopeIndex[2] = root->rope[2]->index;
		root->rope[3] == NULL ? nodeg->RopeIndex[3] = -1 : nodeg->RopeIndex[3] = root->rope[3]->index;
		root->rope[4] == NULL ? nodeg->RopeIndex[4] = -1 : nodeg->RopeIndex[4] = root->rope[4]->index;
		root->rope[5] == NULL ? nodeg->RopeIndex[5] = -1 : nodeg->RopeIndex[5] = root->rope[5]->index;

		nodeg->begin = root->begin;
		nodeg->end = root->end;
		nodeg->box = root->box;
		nodeg->index = root->index;
		nodeg->IsEmpty = root->IsEmpty;
		nodeg->IsLeaf = root->IsLeaf;
		nodeg->PrimCount = root->PrimCount;
		nodeg->SplitPos = root->SplitPos;
		nodeg->Split_Axis = root->Split_Axis;
		_bianli_sz2(root->LeftChild, &(nodeg - nodeg->index)[nodeg->LeftIndex]);
		_bianli_sz2(root->RightChild, &(nodeg - nodeg->index)[nodeg->RightIndex]);
	}
}

KD_Node_V* tree2vector(KD_Node *nodec, int *length)
{
	int index = 0;
	_give_tree_index(nodec, &index);

	KD_Node_V* rootgpu = (KD_Node_V*)malloc(sizeof(KD_Node_V)*index);
	
	_bianli_sz2(nodec, rootgpu);
	*length = index;
	return rootgpu;
}