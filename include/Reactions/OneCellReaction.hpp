#ifndef ONE_CELL_REACTION_HPP
#define ONE_CELL_REACTION_HPP

#include <sundials/sundials_types.h>

#include <limits>

#include "EigenDataTypes.hpp"
#include "Reactions/ReactionSystem.hpp"
#include "Solver.hpp"

std::pair<RowVector, realtype> oneCellReaction(const ReactionSystem& reactionSystem,
                                               RowVector solution,
                                               realtype t_duration,
                                               SolverType solverType,
                                               realtype timeout_seconds = std::numeric_limits<realtype>::infinity());

#endif  // ONE_CELL_REACTION_HPP