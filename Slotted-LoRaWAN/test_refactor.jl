using Test
include("src/SlottedLoRaWAN.jl")
using .SlottedLoRaWAN

@testset "SlottedLoRaWAN Refactor Verification" begin
    @testset "Parameters" begin
        params = create_integrated_params()
        @test params isa IntegratedParameters
        @test params.num_terminals == 400
    end

    @testset "Utilities" begin
        # Dummy data for testing
        time_rx = collect(0.0:0.1:10.0)
        rx_power = ones(length(time_rx)) * 1e-14 # -110 dBm noise floor
        rx_power[40:60] .= 1e-6 # Signal burst
        
        nf = estimate_noise_floor_integrated(rx_power, time_rx)
        @test nf < -100 # Should be around noise level
        
        cross = detect_crossings_integrated(rx_power, time_rx, 1e-9)
        @test !isempty(cross[:start_times])
        @test cross[:start_times][1] ≈ 4.0 atol=0.2
    end

    @testset "Simulation Types" begin
        @test @isdefined(TransmissionRecord)
        @test @isdefined(MACEvent)
        @test TX_START isa EventType
    end

    @testset "Synchronization Engine" begin
        @test @isdefined(perform_hifi_synchronization)
    end

    @testset "MAC Logic" begin
        @test @isdefined(schedule_tx_start!)
    end
end

println("\n✅ All refactor tests passed!")
