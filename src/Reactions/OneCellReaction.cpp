#include "Reactions/OneCellReaction.hpp"

#include <memory>
#include <utility>

#include "EigenDataTypes.hpp"
#include "Observers/SnapshotObserver.hpp"
#include "Process.hpp"
#include "Solver.hpp"
#include "UnitOperations/Volume.hpp"

std::pair<RowVector, realtype> oneCellReaction(const ReactionSystem& reactionSystem,
                                               RowVector solution,
                                               realtype t_duration,
                                               SolverType solverType,
                                               realtype timeout_seconds) {
    auto volume = std::make_shared<Volume>(reactionSystem,
                  1.0);  // Volume of 1 m³ (solution is concentration, y contains masses/amounts -> no conversion necessary)
    volume->y = solution;

    Process process{reactionSystem.componentSystem, {volume}};
    Solver solver{process, solverType};

    auto obs = std::make_shared<SnapshotObserver>(t_duration, volume->all(), true);
    solver.add(obs);

    solver.solve(t_duration, timeout_seconds);

    assert(obs->snapshot.rows() == 1 && "Volume must not have more than 1 cell!");
    RowVector reacted_solution = obs->snapshot;

    auto error = obs->error;
    assert(error != -1.0 && "Error was not computed!");

    return std::make_pair(reacted_solution, error);
}
