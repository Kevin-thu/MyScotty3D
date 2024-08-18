
#include "bvh.h"
#include "aggregate.h"
#include "instance.h"
#include "tri_mesh.h"

#include <stack>

namespace PT {

template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	//A3T3 - build a bvh

	// Keep these
    nodes.clear();
    primitives = std::move(prims);

    // Construct a BVH from the given vector of primitives and maximum leaf
    // size configuration.

	//TODO
	BBox root_box;
	for (auto& prim: primitives) {
		root_box.enclose(prim.bbox());
	}
	root_idx = new_node(root_box, 0, n_primitives());
	build_node(root_idx, max_leaf_size);
	// for (auto& node : nodes) {
	// 	std::cout << node.start << " " << node.size << " " << node.l << " " << node.r << " [" << node.bbox.min << ", " << node.bbox.max << "]" << std::endl;
	// }
}

template<typename Primitive>
void BVH<Primitive>::build_node(size_t node_id, size_t max_leaf_size) {
	Node& node = nodes[node_id];
	if (node.size <= max_leaf_size) {
		node.l = node.r = node_id;
		return;
	}
	size_t num = (node.size > n_prim_per_bucket) ? n_prim_per_bucket : 1;
	size_t num_buckets = ceil(float(node.size) / num);
	SAHBucketData buckets[num_buckets];
	float best_sah_cost = FLT_MAX;
	size_t best_axis = 0;
	size_t best_partition = 0;
	BBox best_left_bbox, best_right_bbox;
	for (size_t axis = 0; axis < 3; ++axis) {
		// std::cout << "hello" << axis << std::endl;
		std::sort(primitives.begin() + node.start, 
			primitives.begin() + node.start + node.size, 
			[axis](const Primitive& p1, const Primitive& p2) {
				return p1.bbox().center()[axis] < p2.bbox().center()[axis];
			}
		);
		for (size_t i = 0; i < num_buckets; ++i) buckets[i].clear();
		for (size_t i = 0; i < node.size; ++i) {
			buckets[i / num].bbox.enclose((primitives.begin() + node.start + i)->bbox());
			buckets[i / num].num_prims += 1;
		}
		for (size_t i = 1; i < num_buckets; ++i) {
			// std::cout << "Node" << node_id << " Axis" << axis << " Bin" << i << ":" << buckets[i].num_prims << std::endl;
			SAHBucketData left, right;
			for (size_t j = 0; j < i; ++j) {
				left.bbox.enclose(buckets[j].bbox);
				left.num_prims += buckets[j].num_prims;
			}
			for (size_t j = i; j < num_buckets; ++j) {
				right.bbox.enclose(buckets[j].bbox);
				right.num_prims += buckets[j].num_prims;
			}
			float sah_cost = left.num_prims * left.bbox.surface_area() + right.num_prims * right.bbox.surface_area();
			if (sah_cost < best_sah_cost) {
				best_sah_cost = sah_cost;
				best_axis = axis;
				best_partition = i;
				best_left_bbox = left.bbox;
				best_right_bbox = right.bbox;
			}
		}
	}
	std::sort(primitives.begin() + node.start, 
		primitives.begin() + node.start + node.size, 
		[best_axis](const Primitive& p1, const Primitive& p2) {
			return p1.bbox().center()[best_axis] < p2.bbox().center()[best_axis];
		}
	);
	size_t left_size = best_partition * num;
	size_t l = new_node(best_left_bbox, node.start, left_size);
	size_t r = new_node(best_right_bbox, node.start + left_size, node.size - left_size);
	//! 当 new_node 修改 vector 结构后原来的 node 指针/引用会失效！
	nodes[node_id].l = l;
	nodes[node_id].r = r;
	build_node(l, max_leaf_size);
	build_node(r, max_leaf_size);
}

template<typename Primitive> 
Trace BVH<Primitive>::hit(const Ray& ray) const {
	//A3T3 - traverse your BVH

    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

	//TODO: replace this code with a more efficient traversal:
	Trace ret;
	if (!nodes.empty())
		find_closest_hit(ray, &nodes[root_idx], ret);
	return ret;
}

template<typename Primitive> 
void BVH<Primitive>::find_closest_hit(const Ray& ray, const Node* node, Trace& closest) const {
	Vec2 bounds = ray.dist_bounds;
	if (!node->bbox.hit(ray, bounds) || bounds.x > closest.distance) {
		return;
	}
	// std::cout << "root:" << node->bbox.min << node->bbox.max << " " << bounds << " " << closest.distance <<std::endl;
	if (node->is_leaf()) {
		for (size_t prim_id = node->start; prim_id < node->start + node->size; ++prim_id) {
			Trace hit = primitives[prim_id].hit(ray);
			closest = Trace::min(closest, hit);
			// std::cout << "leaf:" << primitives[prim_id].bbox().center() << hit.hit << " " << hit.distance << " " << closest.distance << std::endl;
		}
	}
	else {
		const Node* left = &nodes[node->l];
		const Node* right = &nodes[node->r];
		Vec2 left_bounds = bounds, right_bounds = bounds;
		left->bbox.hit(ray, left_bounds);
		right->bbox.hit(ray, right_bounds);
		// std::cout << "left:" << left->bbox.min << left->bbox.max << " " << left_bounds << " " << closest.distance << std::endl;
		// std::cout << "right:" << right->bbox.min << right->bbox.max << " " << right_bounds << " " << closest.distance << std::endl;
		const Node* first, * second;
		float second_bound;
		if (left_bounds.x <= right_bounds.x) {
			first = left;
			second = right;
			second_bound = right_bounds.x;
		} else {
			first = right;
			second = left;
			second_bound = left_bounds.x;
		}
		find_closest_hit(ray, first, closest);
		if (second_bound < closest.distance) {
			find_closest_hit(ray, second, closest);
		}
	}
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
	build(std::move(prims), max_leaf_size);
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
	nodes.clear();
	return std::move(primitives);
}

template<typename Primitive>
template<typename P>
typename std::enable_if<std::is_copy_assignable_v<P>, BVH<P>>::type BVH<Primitive>::copy() const {
	BVH<Primitive> ret;
	ret.nodes = nodes;
	ret.primitives = primitives;
	ret.root_idx = root_idx;
	return ret;
}

template<typename Primitive> Vec3 BVH<Primitive>::sample(RNG &rng, Vec3 from) const {
	if (primitives.empty()) return {};
	int32_t n = rng.integer(0, static_cast<int32_t>(primitives.size()));
	return primitives[n].sample(rng, from);
}

template<typename Primitive>
float BVH<Primitive>::pdf(Ray ray, const Mat4& T, const Mat4& iT) const {
	if (primitives.empty()) return 0.0f;
	float ret = 0.0f;
	for (auto& prim : primitives) ret += prim.pdf(ray, T, iT);
	return ret / primitives.size();
}

template<typename Primitive> void BVH<Primitive>::clear() {
	nodes.clear();
	primitives.clear();
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
	// A node is a leaf if l == r, since all interior nodes must have distinct children
	return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
	Node n;
	n.bbox = box;
	n.start = start;
	n.size = size;
	n.l = l;
	n.r = r;
	nodes.push_back(n);
	return nodes.size() - 1;
}
 
template<typename Primitive> BBox BVH<Primitive>::bbox() const {
	if(nodes.empty()) return BBox{Vec3{0.0f}, Vec3{0.0f}};
	return nodes[root_idx].bbox;
}

template<typename Primitive> size_t BVH<Primitive>::n_primitives() const {
	return primitives.size();
}

template<typename Primitive>
uint32_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, uint32_t level,
                                   const Mat4& trans) const {

	std::stack<std::pair<size_t, uint32_t>> tstack;
	tstack.push({root_idx, 0u});
	uint32_t max_level = 0u;

	if (nodes.empty()) return max_level;

	while (!tstack.empty()) {

		auto [idx, lvl] = tstack.top();
		max_level = std::max(max_level, lvl);
		const Node& node = nodes[idx];
		tstack.pop();

		Spectrum color = lvl == level ? Spectrum(1.0f, 0.0f, 0.0f) : Spectrum(1.0f);
		GL::Lines& add = lvl == level ? active : lines;

		BBox box = node.bbox;
		box.transform(trans);
		Vec3 min = box.min, max = box.max;

		auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

		edge(min, Vec3{max.x, min.y, min.z});
		edge(min, Vec3{min.x, max.y, min.z});
		edge(min, Vec3{min.x, min.y, max.z});
		edge(max, Vec3{min.x, max.y, max.z});
		edge(max, Vec3{max.x, min.y, max.z});
		edge(max, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
		edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
		edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

		if (!node.is_leaf()) {
			tstack.push({node.l, lvl + 1});
			tstack.push({node.r, lvl + 1});
		} else {
			for (size_t i = node.start; i < node.start + node.size; i++) {
				uint32_t c = primitives[i].visualize(lines, active, level - lvl, trans);
				max_level = std::max(c + lvl, max_level);
			}
		}
	}
	return max_level;
}

template class BVH<Triangle>;
template class BVH<Instance>;
template class BVH<Aggregate>;
template BVH<Triangle> BVH<Triangle>::copy<Triangle>() const;

} // namespace PT
