#include <string>

#include "Logger.hpp"
#include "Reactions/ReactionSystem.hpp"

void ReactionSystem::add(const Reaction& reaction) {
	reactions.push_back(reaction);
	LOG("reaction_system.log", "Added reaction No. " + std::to_string(reactions.size()) + "\n");
}
