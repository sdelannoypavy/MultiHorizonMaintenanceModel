using Plots
using Flux
using InferOpt
using ProgressMeter
using Random
using Revise
using Statistics: mean
using MirrorDescentMultistage
using UnicodePlots
using Distributions

function main()

    seed = 2

    Random.seed!(seed)
    rng = MersenneTwister(seed)

    # Hyperparameters
    lr_start = 0.0002
    nb_samples = 20 # number of perturbed samples
    ε = 0.01  # scale of the perturbation
    α = 0.3 # damping parameter
    Nb_epochs = 5 # inner iterations of SGD
    Nb_it = 50 # outer iterations of mirror descent 
    training_data_nb = 20
    T = 3
    Tmax = 30
    γ_lex = 1.0
    κ_val = 1.0
    nb_simu_gap = 100


    init_state = [0.0, 0.0]
    failure_prob1 = 0.25
    failure_prob2 = 0.25

    loss_type = "singlestage" #available: "multistage", "singlestage"

    model_type = "dense" #available models for multistage: "rnn", "features"
    #avilable models for single stage: "dense"

    if (loss_type == "multistage")&&(model_type == "dense")
        throw(ErrorException("model type model_type$(model_type) and loss type loss_type$(loss_type) are not compatible"))
    end

    if (loss_type == "singlestage")&&(model_type != "dense")
        throw(ErrorException("model type model_type$(model_type) and loss type loss_type$(loss_type) are not compatible"))
    end



    train_data = build_dataset(rng;failure_prob1 = failure_prob1,failure_prob2 = failure_prob2,nb_data = training_data_nb,T=T,initial_state = false)
    #println(train_data)
    #train_data =  [([0.5, 0.5], [0 1; 0 0]),
    #([0.5, 0.5], [0 0; 0 1]),
    #([0.5, 0.5], [0 0; 0 0]),
    #([0.5, 0.5], [0 0; 0 0])]
    #train_data =  [([0.0, 1.0], [0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0])] #tester plusieurs dataset, voir si les comportements sont normaux
    #train_data =  [([0.0, 1.0], [0 0 0 1 0 0 1 0 0; 0 0 1 0 0 0 0 0 0]),
    #([1.0, 1.0], [0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0])] 

    # model corresponds to integer updates in the paper, imitation model to the fractional updates (+ 1/2)
    # the imitation model comes from an imitation learning problem on the primal solutions, thus its name
    #model = Chain(
    #    Dense(nb_features => 10, tanh; bias=true, init=zeros))

    #model = Chain(
    #    Dense(nb_features => 1, tanh; bias=true, init=zeros), X -> dropdims(X; dims=1)
    #)
    

    #5 is dimension of hidden state, can be changed

    if model_type == "rnn"

        model = Chain(
            #Dense(5 => 5, tanh; bias=true, init=zeros),
            #RNN(5 => 3, tanh; bias=true),     
            #Dense(3 => 3, tanh; bias=true, init=zeros),      
            RNN(5 => 32, tanh),     
            Dense(32 => 16, relu), # better results than dense only
            Dense(16 => 3)   
            #RNN(5 => 3, tanh),
            #Dense(3 => 3)   
        )

        model = f64(model)

        #set w to 0

        #for param in trainables(model)[1:3]
        #    param .= 0.0 # if Dense do not initialise Dense parameters to 0, otherwise RNN layer remains hidden during training
        #end

        #for param in trainables(model)[4:4]
        #    param .= rand(rng) # if Dense do not initialise Dense parameters to 0, otherwise RNN layer remains hidden during training
        #end

        #for param in trainables(model)[5:5]
        #    param .= rand(rng) # if Dense do not initialise Dense parameters to 0, otherwise RNN layer remains hidden during training
        #end

    elseif model_type == "dense"
        model = Chain(
            Dense(2 => 3, tanh; bias=true, init=zeros),             
        )
    elseif model_type == "features"
        model = Chain(
            Dense(5 => 3, tanh; bias=true, init=zeros),          
        )
    else
        throw(ErrorException("model_type$(model_type) is not an acceptable model type"))
    end

    
    model_imitation = deepcopy(model)
    model_imitation_anticipative = deepcopy(model)
 
    # fix_horizon because we evaluate on [Tmax] so we want to best policy on [Tmax]
    optim_value = cost(rng, init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap; discount_val = 1.0, fix_horizon = true)

    gap = 100 * (
            average_cost(rng,model, model_type, init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap) - optim_value
        ) / optim_value 

    gaps = Float64[]
    push!(gaps,gap)
    gaps_imit = Float64[]
    push!(gaps_imit, gap)
    gaps_anticipative = Float64[]
    push!(gaps_anticipative, gap)

    augmented_dataset = []

    #println(train_data)

    anticipative_dataset = build_anticipative_dataset(
            train_data,
            T,
            loss_type
    )

    for n_it in 1:Nb_it
        # augment data sets per scenario given the current model and penalization scale
        # it corresponds to the decomposition update : q^{t+1} = \argmin_q S(w^{t}, q)
        augmented_dataset = one_point_per_scenario(
            train_data,
            model_type,
            loss_type;
            model=model,
            ε=ε,
            nb_samples=nb_samples,
            T=T,
            γ = γ_lex,
            κ = κ_val,
            seed = seed
        )


        # define training
        # it corresponds to the coordination update w^{t+1/2} = \argmin_w S(w, q^{t+1})
        model_imitation, _ = train_model(
            model_imitation,
            augmented_dataset,
            model_type;
            nb_epochs=Nb_epochs,
            ε=ε,
            nb_samples=nb_samples,
            seed=seed,
            lr_start,
            T=T,
            γ=γ_lex
        )

        model_imitation_anticipative, _ = train_model(
            model_imitation_anticipative,
            anticipative_dataset,
            model_type;
            nb_epochs=Nb_epochs,
            ε=ε,
            nb_samples=nb_samples,
            seed=seed,
            lr_start,
            T=T,
            γ=γ_lex
        )
        # damped update 
        # it corresponds to w^{t+1} = α w^{t+1/2} + (1-α) w^{t}

        nb_pramams = length(trainables(model)) 

        #println(trainables(model))

        for p in 1:nb_pramams
            trainables(model)[p].+= 
                α * (trainables(model_imitation)[p] - trainables(model)[p])
        end

        #println(trainables(model))

        Flux.reset!(model)
        # gap = eval_gap(model; init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap)

        gap = 100 * (
            average_cost(rng, model, model_type, init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap) - optim_value
        ) / optim_value 

        gap_imit = 100 * (
            average_cost(rng, model_imitation, model_type, init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap) - optim_value
        ) / optim_value 

        gap_anticipative = 100 * (
            average_cost(rng, model_imitation_anticipative, model_type, init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap) - optim_value
        ) / optim_value 
        
        push!(gaps, gap)
        push!(gaps_imit, gap_imit)
        push!(gaps_anticipative, gap_anticipative)


    end

    if (Nb_it > 1)
        xticks_positions = range(0, stop=Nb_it, length=Nb_it)
        xticks_positions = 0:Nb_it

        plot = Plots.plot(0:Nb_it, gaps,
            label="Gaps", 
            xlabel="Primal-dual iteration t", 
            ylabel="Average gap",
            title=" ",
            color=:black,  # Couleur noire pour la ligne,
            xticks=xticks_positions,
        )
        
        # Sauvegarder le graphique en PDF
        mkpath("scripts/results")
        Plots.savefig(plot, joinpath(
            "scripts",
            "results",
            "mirror_descent_ϵ$(ε)_α$(α)_nsa$(nb_samples)_Nit$(Nb_it)_Nepochs$(Nb_epochs)_lrstart$(lr_start)_κ$(κ_val)_failure_prob1$(failure_prob1)_failure_prob2$(failure_prob2)_γ_lex$(γ_lex)_T$(T)_Tmax$(Tmax).pdf",
        ))
    end

    return model, gaps, gaps_imit, gaps_anticipative, augmented_dataset

end

# to display policies: multistage_maximizer(vec(model(input)))[:,1]


function main_count_right()

    nb_test = 50

    Nb_it = 15

    res = [0.0 for _ in 1:Nb_it]

    for i in 1:nb_test

        seed = i

        Random.seed!(seed)
        rng = MersenneTwister(seed)

        # Hyperparameters
        lr_start = 0.01
        nb_samples = 20 # number of perturbed samples
        ε = 0.01  # scale of the perturbation
        α = 0.3 # damping parameter
        Nb_epochs = 50 # inner iterations of SGD
        T = 3
        γ_lex = 1.0
        κ_val = 1.0


        init_state = [0.0, 0.0]
        failure_prob1 = 0.3
        failure_prob2 = 0.3

        loss_type = "multistage" #available: "multistage", "singlestage"

        model_type = "features" #available models for multistage: "rnn", "features"
        #avilable models for single stage: "dense"

        if (loss_type == "multistage")&&(model_type == "dense")
            throw(ErrorException("model type model_type$(model_type) and loss type loss_type$(loss_type) are not compatible"))
        end

        if (loss_type == "singlestage")&&(model_type != "dense")
            throw(ErrorException("model type model_type$(model_type) and loss type loss_type$(loss_type) are not compatible"))
        end



        #train_data = build_dataset(rng;failure_prob1 = failure_prob1,failure_prob2 = failure_prob2,nb_data = training_data_nb,T=T,initial_state = false)
        #println(train_data)
        train_data =  [([0.5, 0.5], [0 1; 0 0]),
        ([0.5, 0.5], [0 0; 0 1]),
        ([0.5, 0.5], [0 0; 0 0])]
        #train_data =  [([0.0, 1.0], [0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 0 0])] #tester plusieurs dataset, voir si les comportements sont normaux
        #train_data =  [([0.0, 1.0], [0 0 0 1 0 0 1 0 0; 0 0 1 0 0 0 0 0 0]),
        #([1.0, 1.0], [0 0 0 0 0 0 0 0 0; 0 0 0 0 0 0 0 1 0])] 

        # model corresponds to integer updates in the paper, imitation model to the fractional updates (+ 1/2)
        # the imitation model comes from an imitation learning problem on the primal solutions, thus its name
        #model = Chain(
        #    Dense(nb_features => 10, tanh; bias=true, init=zeros))

        #model = Chain(
        #    Dense(nb_features => 1, tanh; bias=true, init=zeros), X -> dropdims(X; dims=1)
        #)
        

        #5 is dimension of hidden state, can be changed

        if model_type == "rnn"

            model = Chain(
                #Dense(5 => 5, tanh; bias=true, init=zeros),
                #RNN(5 => 3, tanh; bias=true),     
                #Dense(3 => 3, tanh; bias=true, init=zeros),      
                RNN(5 => 32, tanh),     
                Dense(32 => 16, relu), # better results than dense only
                Dense(16 => 3)   
                #RNN(5 => 3, tanh),
                #Dense(3 => 3)   
            )

            model = f64(model)

            #set w to 0

            #for param in trainables(model)[1:3]
            #    param .= 0.0 # if Dense do not initialise Dense parameters to 0, otherwise RNN layer remains hidden during training
            #end

            #for param in trainables(model)[4:4]
            #    param .= rand(rng) # if Dense do not initialise Dense parameters to 0, otherwise RNN layer remains hidden during training
            #end

            #for param in trainables(model)[5:5]
            #    param .= rand(rng) # if Dense do not initialise Dense parameters to 0, otherwise RNN layer remains hidden during training
            #end

        elseif model_type == "dense"
            model = Chain(
                Dense(2 => 3, tanh; bias=true, init=zeros),             
            )
        elseif model_type == "features"
            model = Chain(
                Dense(5 => 3, tanh; bias=true, init=zeros),          
            )
        else
            throw(ErrorException("model_type$(model_type) is not an acceptable model type"))
        end

        
        model_imitation = deepcopy(model)

        augmented_dataset = []

        #println(train_data)

        for n_it in 1:Nb_it
            # augment data sets per scenario given the current model and penalization scale
            # it corresponds to the decomposition update : q^{t+1} = \argmin_q S(w^{t}, q)
            augmented_dataset = one_point_per_scenario(
                train_data,
                model_type,
                loss_type;
                model=model,
                ε=ε,
                nb_samples=nb_samples,
                T=T,
                γ = γ_lex,
                κ = κ_val,
                seed = seed
            )
            # define training
            # it corresponds to the coordination update w^{t+1/2} = \argmin_w S(w, q^{t+1})
            model_imitation, _ = train_model(
                model_imitation,
                augmented_dataset,
                model_type;
                nb_epochs=Nb_epochs,
                ε=ε,
                nb_samples=nb_samples,
                seed=seed,
                lr_start,
                T=T,
                γ=γ_lex
            )
            # damped update 
            # it corresponds to w^{t+1} = α w^{t+1/2} + (1-α) w^{t}

            nb_pramams = length(trainables(model)) 

            #println(trainables(model))

            for p in 1:nb_pramams
                trainables(model)[p].+= 
                    α * (trainables(model_imitation)[p] - trainables(model)[p])
            end

            #println(trainables(model))

            Flux.reset!(model)
            # gap = eval_gap(model; init_state, failure_prob1, failure_prob2, Tmax, nb_simu_gap)

            if (multistage_maximizer(vec(model([0.5,0.5,1,0,0])))[:,1] != [0.0, 0.0, 1.0])
                res[n_it] += 1
            end

        end
    
    end

    return res ./ nb_test

end