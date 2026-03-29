module SlottedLoRaWAN

# Order matters for dependencies between modules
# 1. Types and Basic Logic
include("modules/parameters.jl")
include("modules/utilities.jl")
include("modules/simulation_types.jl")
include("modules/lora_airtime.jl")
include("modules/path_loss.jl")
include("modules/shadowing.jl")
include("modules/noise_generation.jl")
include("modules/signal_generation.jl")
include("modules/local_clock.jl")
include("modules/collision_detection.jl")
include("modules/terminal_deployment.jl")
include("modules/packet_generation.jl")

# 2. Advanced Simulation Logic (using above modules)
include("modules/synchronization_engine.jl")
include("modules/mac_logic.jl")

# Re-export key components
using .Parameters, .Utilities, .SimulationTypes, .LoraAirtime
using .PathLoss, .Shadowing, .NoiseGeneration, .SignalGeneration
using .LocalClock, .CollisionDetection, .TerminalDeployment, .PacketGeneration
using .SynchronizationEngine, .MACLogic

# Forward exports
export Parameters, Utilities, SimulationTypes, LoraAirtime
export PathLoss, Shadowing, NoiseGeneration, SignalGeneration
export LocalClock, CollisionDetection, TerminalDeployment, PacketGeneration
export SynchronizationEngine, MACLogic

# Flatten most important common types and functions for easier use
export IntegratedParameters, create_integrated_params
export estimate_noise_floor_integrated, detect_crossings_integrated, debounce_integrated
export SyncResult, CandidateSlot, TransmissionRecord, MACEvent, EventType, TX_START, TX_END, ACK_CHECK, PACKET_GEN
export perform_hifi_synchronization
export schedule_tx_start!
export deploy_terminals, calculate_lora_airtime, create_lora_params
export SignalParameters, calculate_beacon_times
export PathLossParameters, calculate_path_loss, calculate_path_loss_simple
export ShadowingParameters, calculate_shadowing, calculate_correlated_shadowing
export TerminalDeploymentParameters, TerminalInfo

end # module
