
#include "transform.h"

Mat4 Transform::local_to_parent() const {
	return Mat4::translate(translation) * rotation.to_mat() * Mat4::scale(scale);
}

Mat4 Transform::parent_to_local() const {
	return Mat4::scale(1.0f / scale) * rotation.inverse().to_mat() * Mat4::translate(-translation);
}

Mat4 Transform::local_to_world() const {
	// A1T1: local_to_world
	//don't use Mat4::inverse() in your code.
	Mat4 r;
	const Transform* cur = this;
	while (std::shared_ptr <Transform> parent_ = cur->parent.lock()) {
		r = cur->local_to_parent() * r;
		cur = parent_.get();
	}
	return cur->local_to_parent() * r;
}

Mat4 Transform::world_to_local() const {
	// A1T1: world_to_local
	//don't use Mat4::inverse() in your code.
	Mat4 r;
	const Transform* cur = this;
	while (std::shared_ptr <Transform> parent_ = cur->parent.lock()) {
		r *= cur->parent_to_local();
		cur = parent_.get();
	}
	return r * cur->parent_to_local();
}

bool operator!=(const Transform& a, const Transform& b) {
	return a.parent.lock() != b.parent.lock() || a.translation != b.translation ||
	       a.rotation != b.rotation || a.scale != b.scale;
}
