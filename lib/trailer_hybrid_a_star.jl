#
# Hybrid A* path planning for trailer
#
# author: Atsushi Sakai(@Atsushi_twi)
#

module trailer_hybrid_a_star

using PyPlot
using DataFrames
using NearestNeighbors
using DataStructures 
using Base
using Printf
include("rs_path.jl")
include("grid_a_star.jl")
include("trailerlib.jl")

const XY_GRID_RESOLUTION = 2.0 #[m]
const YAW_GRID_RESOLUTION = deg2rad(15.0) #[rad]
const N_PATHS_NEEDED = 1
const GOAL_TYAW_TH = deg2rad(5.0) #[rad]
const MOTION_RESOLUTION = 0.1 #[m] path interporate resolution
const N_STEER = 20.0 # number of steer command
const EXTEND_AREA= 5.0 #[m] map extend length
const SKIP_COLLISION_CHECK= 4 # skip number for collision check
 
const SB_COST = 100.0 # switch back penalty cost
const BACK_COST = 5.0 # backward penalty cost
const STEER_CHANGE_COST = 5.0 # steer angle change penalty cost
const STEER_COST = 1.0 # steer angle change penalty cost
const JACKKNIF_COST= 200.0 # Jackknif cost
const H_COST = 5.0 # Heuristic cost
 
const WB = trailerlib.WB #[m] Wheel base
const LT = trailerlib.LT #[m] length of trailer
const MAX_STEER = trailerlib.MAX_STEER #[rad] maximum steering angle


struct Node
    xind::Int64 #x index
    yind::Int64 #y index
    yawind::Int64 #yaw index
    direction::Bool # moving direction forword:true, backword:false
    x::Array{Float64} # x position [m]
    y::Array{Float64} # y position [m]
    yaw::Array{Float64} # yaw angle [rad]
    yaw1::Array{Float64} # trailer yaw angle [rad]
    directions::Array{Bool} # directions of each points forward: true, backward:false
    steer::Float64 # steer input
    cost::Float64 # cost
    pind::Int64 # parent index
end

function Base.show(io::IO, node::Node)
    println(io, "Node( xind: ", node.xind, "  yind: ", node.yind
                    , "  yawind: ", node.yawind, "  direction: ", node.direction ? "forward" : "backward" 
                    , "  steer: ", node.steer, "  cost: ", node.cost, "  pind: ", node.pind, "poses size: ", length(node.y), ")")
    # println(io, )
    # println(io)"
    # println(io)
    # println(io)
    # # println(io, "  x: ", node.x)
    # # println(io, "  y: ", node.y)
    # # println(io, "  yaw: ", node.yaw)
    # # println(io, "  yaw1: ", node.yaw1)
    # # println(io, "  directions: ", node.directions)
    # println(io, "  steer: ", node.steer)
    # println(io, "  cost: ", node.cost)
    # println(io, "  pind: ", node.pind)
    # print(io, ")")
end

struct Config # config struct for hybrid A* DB
    minx::Int64
    miny::Int64
    minyaw::Int64
    minyawt::Int64
    maxx::Int64
    maxy::Int64
    maxyaw::Int64
    maxyawt::Int64
    xw::Int64
    yw::Int64
    yaww::Int64
    yawtw::Int64
    xyreso::Float64
    yawreso::Float64
end

mutable struct Path
    x::Array{Float64} # x position [m]
    y::Array{Float64} # y position [m]
    yaw::Array{Float64} # yaw angle [rad]
    yaw1::Array{Float64} # trailer angle [rad]
    direction::Array{Bool} # direction forward: true, back false
    cost::Float64 # cost
end

Base.length(p::Path) = length(p.x)

function calc_hybrid_astar_path(sx::Float64, sy::Float64, syaw::Float64, syaw1::Float64,
                                gx::Float64, gy::Float64, gyaw::Float64, gyaw1::Float64,
                                ox::Array{Float64}, oy::Array{Float64},
                                xyreso::Float64, yawreso::Float64,
                                )
    """
    sx: start x position [m]
    sy: start y position [m]
    gx: goal x position [m]
    gx: goal x position [m]
    ox: x position list of Obstacles [m]
    oy: y position list of Obstacles [m]
    xyreso: grid resolution [m]
    yawreso: yaw angle resolution [rad]
    """

    syaw, gyaw = rs_path.pi_2_pi(syaw), rs_path.pi_2_pi(gyaw)
    ox, oy = ox[:], oy[:]

    kdtree = KDTree(hcat(ox, oy)')
    # println(" ==========s=== ", sx, sy, syaw, syaw1, xyreso, yawreso);
    # println(" ==========g=== ", gx, gy, gyaw, gyaw1, xyreso, yawreso);
            
    c = calc_config(ox, oy, xyreso, yawreso)
    nstart = Node(round(Int64,sx/xyreso), round(Int64,sy/xyreso), round(Int64, syaw/yawreso),true,[sx],[sy],[syaw],[syaw1],[true],0.0,0.0, -1)
    ngoal = Node(round(Int64,gx/xyreso), round(Int64,gy/xyreso), round(Int64,gyaw/yawreso),true,[gx],[gy],[gyaw],[gyaw1],[true],0.0,0.0, -1)
    
    # println("nstart: ", nstart);
    # println("ngoal: ", nstart);
    h_dp = calc_holonomic_with_obstacle_heuristic(ngoal, ox, oy, xyreso)

    open, closed = Dict{Int64, Node}(), Dict{Int64, Node}()
    fnode = nothing
    open[calc_index(nstart, c)] = nstart
    pq = PriorityQueue{Int64,Float64}()
    enqueue!(pq, calc_index(nstart, c), calc_cost(nstart, h_dp, ngoal, c))

    u, d = calc_motion_inputs()
    nmotion = length(u)
    # println("u d ", nmotion,  " " , length(d));
    # println(u);
    # println(d);
    counter = 0 
    while true
        # break;
        if length(open) == 0
            println("Error: Cannot find path, No open set")
            return []
        end

        c_id = dequeue!(pq)
        current = open[c_id]

        #move current node from open to closed
        delete!(open, c_id)
        closed[c_id] = current
        # println(" c_id ", c_id)
        # println("------------------------node info------------------------")
        # println(" current ", current)
        # println(" ngoal ", ngoal)
        # println(" gyaw1 ", gyaw1)
        isupdated, fpath = update_node_with_analystic_expantion(current, ngoal, c, ox, oy, kdtree, gyaw1)
        # println(" isupdated ", isupdated)
        # println(" fpath ", fpath)
        # break
        if isupdated # found
            fnode = fpath
            break
        end
        # println(" -------current poses----------- ", current.x, " ", current.y, " ", current.yaw, " ", current.yaw1)
        inityaw1 = current.yaw1[1]
        # println("-----------------------------------")
        for i in 1:nmotion
            node = calc_next_node(current, c_id, u[i], d[i], c)
           
            if !verify_index(node, c, ox, oy, inityaw1, kdtree) continue end

            node_ind = calc_index(node, c)
            # println("-------------node info----------------------")
            # println(" node_ind ", node_ind)
            # println(" node ", node)
            # println("-------------end---------------------")
            # If it is already in the closed set, skip it
            if haskey(closed, node_ind)  continue end

            if !haskey(open, node_ind)
                open[node_ind] = node
                cost_node = calc_cost(node, h_dp, ngoal, c)
                # println(" node_ind ", node_ind, " cost_node ", cost_node)
                # break
                enqueue!(pq, node_ind, cost_node)
            else
                if open[node_ind].cost > node.cost
                    # If so, update the node to have a new parent
                    # println(" open ", open[node_ind].cost, " cost ", node.cost, " id: ", node.pind)
                    open[node_ind] = node
                end
            end

            # counter = counter + 1
            # if(counter == 5)
            #     break
            # end
        end
        # println("-----------------------------------")
        # break
        # counter = counter + 1
        # if(counter == 2)
        #     break
        # end
    end

    println("final expand node:", length(open) + length(closed))

    path = get_final_path(closed, fnode, nstart, c)

    return path
end


function update_node_with_analystic_expantion(current::Node,
                                             ngoal::Node,
                                             c::Config,
                                             ox,
                                             oy,
                                             kdtree,
                                             gyaw1::Float64
                                            )

    apath = analystic_expantion(current, ngoal, c, ox, oy, kdtree)
    # println(" apath ", apath)
    if apath != nothing
        # println("----analystic_expantion---- ")
        fx = apath.x[2:end]
        fy = apath.y[2:end]
        fyaw =  apath.yaw[2:end]
        steps = MOTION_RESOLUTION*apath.directions
        # println("steps: ", steps)
        # println(" yaw: ", apath.yaw)
        yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(apath.x, apath.y, apath.yaw, current.yaw1[end], steps)
        # println(" yaw1 ", yaw1[end], " gyaw1 ", gyaw1, " GOAL_TYAW_TH ", GOAL_TYAW_TH, " pi_2_pi ",  rs_path.pi_2_pi(yaw1[end] - gyaw1), " yaw1 ", current.yaw1[end])
        if abs(rs_path.pi_2_pi(yaw1[end] - gyaw1)) >= GOAL_TYAW_TH
            return false, nothing #no update
        end
        fcost = current.cost + calc_rs_path_cost(apath, yaw1)
        fyaw1 = yaw1[2:end]
        fpind = calc_index(current, c)

        fd = Bool[]
        for d in apath.directions[2:end]
            if d >= 0
                push!(fd, true)
            else
                push!(fd, false)
            end
        end

        fsteer = 0.0

        fpath = Node(current.xind, current.yind, current.yawind, current.direction, fx, fy, fyaw, fyaw1, fd, fsteer, fcost, fpind)

        return true, fpath
    end

    return false, nothing #no update
end


function calc_rs_path_cost(rspath::rs_path.Path, yaw1)

    cost = 0.0
    for l in rspath.lengths
        if l >= 0 # forward
            cost += l
        else # back
            cost += abs(l) * BACK_COST
        end
    end

    # println(" 1cost ", cost )

    # swich back penalty
    for i in 1:length(rspath.lengths) - 1
        if rspath.lengths[i] * rspath.lengths[i+1] < 0.0 # switch back
            cost += SB_COST
        end
    end
    # println(" 2cost ", cost )

    # steer penalyty
    for ctype in rspath.ctypes
        if ctype != "S" # curve
            cost += STEER_COST*abs(MAX_STEER)
        end
    end
    # println(" 3cost ", cost )
    # ==steer change penalty
    # calc steer profile
    nctypes = length(rspath.ctypes)
    ulist = fill(0.0, nctypes)
    for i in 1:nctypes
        if rspath.ctypes[i] == "R" 
            ulist[i] = - MAX_STEER
        elseif rspath.ctypes[i] == "L"
            ulist[i] = MAX_STEER
        end
    end
    # println(" 4cost ", cost )
    for i in 1:length(rspath.ctypes) - 1
        cost += STEER_CHANGE_COST*abs(ulist[i+1] - ulist[i])
    end
    # println(" 5cost ", cost )
    # println(" rspath.yaw ", rspath.yaw)
    # println(" yaw1 ", yaw1)
    angle_diff = sum(abs.(rs_path.pi_2_pi.(rspath.yaw-yaw1)))
    cost += JACKKNIF_COST * angle_diff
    # println(" 6cost ", cost, " angle_diff ", angle_diff)
    return cost
end


function analystic_expantion(n::Node, ngoal::Node, c::Config, ox, oy, kdtree)

    sx = n.x[end]
    sy = n.y[end]
    syaw = n.yaw[end]

    max_curvature = tan(MAX_STEER)/WB
    paths = rs_path.calc_paths(sx,sy,syaw,ngoal.x[end], ngoal.y[end], ngoal.yaw[end],
                                   max_curvature, step_size=MOTION_RESOLUTION)
    # println(" 1analystic_expantion 1")
    if length(paths) == 0 
        return nothing
    end

    # println(" 1analystic_expantion paths ", length(paths));

    pathqueue = PriorityQueue{rs_path.Path, Float64}()
    for path in paths
        steps = MOTION_RESOLUTION*path.directions
        yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(path.x, path.y, path.yaw, n.yaw1[end], steps)
        cost = calc_rs_path_cost(path, yaw1)
        enqueue!(pathqueue, path, cost)
        # println(" 1analystic_expantion path.cost ", cost);
    end

    for i in length(pathqueue)
        path = dequeue!(pathqueue)

        steps = MOTION_RESOLUTION*path.directions
        yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(path.x, path.y, path.yaw, n.yaw1[end], steps)
        ind = 1:SKIP_COLLISION_CHECK:length(path.x)
        if trailerlib.check_trailer_collision(ox, oy, path.x[ind], path.y[ind], path.yaw[ind], yaw1[ind], kdtree = kdtree)
            # plot(path.x, path.y, "-^b")
            return path # path is ok
        end
    end
    # println(" 1analystic_expantion 2")

    return nothing
end


function calc_motion_inputs()

    up = [i for i in MAX_STEER/N_STEER:MAX_STEER/N_STEER:MAX_STEER]
    u = vcat([0.0], [i for i in up], [-i for i in up]) 
    d = vcat([1.0 for i in 1:length(u)], [-1.0 for i in 1:length(u)]) 
    u = vcat(u,u)

    return u, d
end


function verify_index(node::Node, c::Config, ox, oy, inityaw1, kdtree)

    # overflow map
    if (node.xind - c.minx) >= c.xw
        return false
    elseif (node.xind - c.minx) <= 0
        return false
    end
    if (node.yind - c.miny) >= c.yw
        return false
    elseif (node.yind - c.miny) <= 0
        return false
    end

    # check collisiton
    steps = MOTION_RESOLUTION*node.directions
    yaw1 = trailerlib.calc_trailer_yaw_from_xyyaw(node.x, node.y, node.yaw, inityaw1, steps)
    ind = 1:SKIP_COLLISION_CHECK:length(node.x)
    if !trailerlib.check_trailer_collision(ox, oy, node.x[ind], node.y[ind], node.yaw[ind], yaw1[ind], kdtree = kdtree)
        return false
    end

    return true #index is ok"
end


function calc_next_node(current::Node, c_id::Int64,
                        u::Float64, d::Float64, 
                        c::Config,
                        )

    arc_l = XY_GRID_RESOLUTION*1.5

    nlist = floor(Int64, arc_l/MOTION_RESOLUTION)+1
    xlist = fill(0.0, nlist)
    ylist = fill(0.0, nlist)
    yawlist = fill(0.0, nlist)
    yaw1list = fill(0.0, nlist)
    # println(" current pose ", current.x, " nlist ", nlist)
    # println(" current pose ", current.y)
    xlist[1] = current.x[end] + d * MOTION_RESOLUTION*cos(current.yaw[end])
    ylist[1] = current.y[end] + d * MOTION_RESOLUTION*sin(current.yaw[end])
    yawlist[1] = rs_path.pi_2_pi(current.yaw[end] + d*MOTION_RESOLUTION/WB * tan(u))
    yaw1list[1] = rs_path.pi_2_pi(current.yaw1[end] + d*MOTION_RESOLUTION/LT*sin(current.yaw[end]-current.yaw1[end]))
 
    for i in 1:(nlist-1)
        xlist[i+1] = xlist[i] + d * MOTION_RESOLUTION*cos(yawlist[i])
        ylist[i+1] = ylist[i] + d * MOTION_RESOLUTION*sin(yawlist[i])
        yawlist[i+1] = rs_path.pi_2_pi(yawlist[i] + d*MOTION_RESOLUTION/WB * tan(u))
        yaw1list[i+1] = rs_path.pi_2_pi(yaw1list[i] + d*MOTION_RESOLUTION/LT*sin(yawlist[i]-yaw1list[i]))
    end

    # println(" xlist ", xlist)
    # println(" ylist ", ylist)
    # println(" yawlist ", yawlist)
    # println(" yaw1list ", yaw1list)
 
    xind = round(Int64, xlist[end]/c.xyreso)
    yind = round(Int64, ylist[end]/c.xyreso)
    yawind = round(Int64, yawlist[end]/c.yawreso)

    # println(" xind ", xind, " yind ", yind, " yawind ", yawind) 

    # println(" =========================cost=========================== ") 
    addedcost = 0.0
    if d > 0
        direction = true
        addedcost += abs(arc_l)
    else
        direction = false
        addedcost += abs(arc_l) * BACK_COST
    end
    # println(" =========== addedcost1 ", addedcost) 


    # swich back penalty
    if direction != current.direction # switch back penalty
        addedcost += SB_COST
    end

    # println(" =========== addedcost2 ", addedcost) 

    # steer penalyty
    addedcost += STEER_COST*abs(u)

    # println(" =========== addedcost3 ", addedcost) 

    # steer change penalty
    addedcost += STEER_CHANGE_COST*abs(current.steer - u)
    # println(" =========== addedcost4 ", addedcost) 

    # println(" =========== yawlist ", yawlist, " ----- ", yaw1list) 
    angle_diff = sum(abs.(rs_path.pi_2_pi.(yawlist-yaw1list)))
    # jacknif cost
    addedcost += JACKKNIF_COST * angle_diff
    # println(" =========== addedcost5 ", addedcost, " angle_diff ", angle_diff) 

    cost = current.cost + addedcost 
    # println(" =========== addedcost6 ", addedcost) 

    directions = [direction for i in 1:length(xlist)]
    node = Node(xind, yind, yawind, direction, xlist, ylist, yawlist, yaw1list, directions, u, cost, c_id)

    return node
end


function is_same_grid(node1::Node,node2::Node)

    if node1.xind != node2.xind
        return false
    end
    if node1.yind != node2.yind
        return false
    end
    if node1.yawind != node2.yawind
        return false
    end

    return true

end


function calc_index(node::Node, c::Config)
    ind = (node.yawind - c.minyaw)*c.xw*c.yw+(node.yind - c.miny)*c.xw + (node.xind - c.minx)

    # 4D grid
    yaw1ind = round(Int64, node.yaw1[end]/c.yawreso)
    ind += (yaw1ind - c.minyawt) *c.xw*c.yw*c.yaww

    if ind <= 0
        println("Error(calc_index):", ind)
    end
    return ind
end


function calc_holonomic_with_obstacle_heuristic(gnode::Node, ox::Array{Float64}, oy::Array{Float64}, xyreso::Float64)
    h_dp = grid_a_star.calc_dist_policy(gnode.x[end], gnode.y[end], ox, oy, xyreso, 1.0)
    return h_dp
end


function calc_config(ox::Array{Float64}, oy::Array{Float64},
                     xyreso::Float64, yawreso::Float64
                     )

    min_x_m = minimum(ox) - EXTEND_AREA
    min_y_m = minimum(oy) - EXTEND_AREA
    max_x_m = maximum(ox) + EXTEND_AREA
    max_y_m = maximum(oy) + EXTEND_AREA

    push!(ox, min_x_m)
    push!(oy, min_y_m)
    push!(ox, max_x_m)
    push!(oy, max_y_m)

    minx = round(Int64, min_x_m/xyreso)
    miny = round(Int64, min_y_m/xyreso)
    maxx = round(Int64, max_x_m/xyreso)
    maxy = round(Int64, max_y_m/xyreso)

    
    xw = round(Int64,(maxx - minx))
    yw = round(Int64,(maxy - miny))
    println("   minx miny maxx maxy xw yw ", minx, " ", miny, " ", maxx , " ", maxy, " ", xw, " ", yw )

    minyaw = round(Int64, - pi/yawreso) - 1
    maxyaw = round(Int64, pi/yawreso)
    yaww = round(Int64,(maxyaw - minyaw))

    minyawt = minyaw
    maxyawt = maxyaw
    yawtw = yaww

    config = Config(minx, miny, minyaw, minyawt, maxx, maxy, maxyaw, maxyawt, xw, yw, yaww, yawtw,
                    xyreso, yawreso)

    return config
end


function get_final_path(closed::Dict{Int64, Node},
                        ngoal::Node,
                        nstart::Node,
                        c::Config)

    rx, ry, ryaw = Array{Float64}(reverse(ngoal.x)),Array{Float64}(reverse(ngoal.y)),Array{Float64}(reverse(ngoal.yaw))
    ryaw1 = Array{Float64}(reverse(ngoal.yaw1))
    direction = Array{Float64}(reverse(ngoal.directions))
    nid = ngoal.pind
    finalcost = ngoal.cost

    println("Final path info ", nid , " finalcost ", finalcost)

    while true
        n = closed[nid]
        rx = vcat(rx, reverse(n.x))
        ry = vcat(ry, reverse(n.y))
        ryaw = vcat(ryaw, reverse(n.yaw))
        ryaw1 = vcat(ryaw1, reverse(n.yaw1))
        direction = vcat(direction, reverse(n.directions))
        nid = n.pind
        if is_same_grid(n, nstart)
            break
        end
    end

    rx = reverse(rx)
    ry = reverse(ry)
    ryaw = reverse(ryaw)
    ryaw1 = reverse(ryaw1)
    direction = reverse(direction)

    # adjuct first direction
    direction[1] = direction[2]

    path = Path(rx, ry, ryaw, ryaw1, direction, finalcost)

    # println(" final_path length ",length(rx));
    # println(" final_path length ", ryaw);

    return path
end


function calc_cost(n::Node, h_dp::Array{Float64}, ngoal::Node, c::Config)
    
    # formatted_A = map(elem -> @sprintf("%.4f", elem), h_dp)
    # # println(h_dp)
    # # Print the formatted matrix
    # for row in eachrow(formatted_A)
    #     println(join(row, " "))
    # end
    # for row in h_dp
    #     for elem in row
    #         println(round.(elem; sigdigits=4))
    #     end
    #     println()  # Newline at the end of each row
    # end
    # println(" n.xind ", n.xind, " n.yind ", n.yind, " ", n.xind - c.minx, " ", n.yind - c.miny)
    if(n.xind - c.minx < 0)
        println("=======================x======================================================")
    end
    if(n.yind - c.miny < 0)
        println("======================y=======================================================")
    end
    total_cost = (n.cost + H_COST*h_dp[n.xind - c.minx, n.yind - c.miny])
    # println(" n.cost ", n.cost, " hp ", h_dp[n.xind - c.minx, n.yind - c.miny], " total: ", total_cost)
   return total_cost

end


function main()
    println(PROGRAM_FILE," start!!")

    sx = 14.0  # [m]
    sy = 10.0  # [m]
    syaw0 = deg2rad(00.0)
    syaw1 = deg2rad(00.0)

    gx = 0.0  # [m]
    gy = 0.0  # [m]
    gyaw0 = deg2rad(90.0)
    gyaw1 = deg2rad(90.0)

    ox = Float64[]
    oy = Float64[]

    for i in -25:25
        push!(ox, Float64(i))
        push!(oy, 15.0)
    end
    for i in -25:-4
        push!(ox, Float64(i))
        push!(oy, 4.0)
    end
    for i in -15:4
        push!(ox, -4.0)
        push!(oy, Float64(i))
    end
    for i in -15:4
        push!(ox, 4.0)
        push!(oy, Float64(i))
    end
    for i in 4:25
        push!(ox, Float64(i))
        push!(oy, 4.0)
    end
    for i in -4:4
        push!(ox, Float64(i))
        push!(oy, -15.0)
    end

    oox = ox[:]
    ooy = oy[:]

    @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    plot(oox, ooy, ".k")
    trailerlib.plot_trailer(sx, sy, syaw0, syaw1, 0.0)
    trailerlib.plot_trailer(gx, gy, gyaw0, gyaw1, 0.0)
    x = path.x
    y = path.y
    yaw = path.yaw
    yaw1 = path.yaw1
    direction = path.direction
    println("direction: ", direction)
    steer = 0.0
    for ii in 1:length(x)
        cla()
        plot(oox, ooy, ".k")
        plot(x, y, "-r", label="Hybrid A* path")

        if ii < length(x)-1
            k = (yaw[ii+1] - yaw[ii])/MOTION_RESOLUTION
            println("k: ", direction[ii], ",", ii)
            if !direction[ii]
                # println("index: ", direction[ii])
                k *= -1
            end
            steer = Base.atan(WB*k, 1.0)
        else
            println("index: ", ii)
            steer = 0.0
        end
        trailerlib.plot_trailer.(x[ii], y[ii], yaw[ii], yaw1[ii], steer)
        grid(true)
        axis("equal")
        pause(0.0001)
        println("steer: ", x[ii], "," ,y[ii], "," ,yaw[ii],  ",",yaw1[ii], "," ,steer)
    end
    println("Done")
    axis("equal")
    show()

    println(PROGRAM_FILE," Done!!")
end


function test()
    # println("Test Start !!!")

    gx = 0.0  # [m]
    gy = 0.0  # [m]
    gyaw0 = deg2rad(90.0)
    gyaw1 = deg2rad(90.0)


    ox = Float64[]
    oy = Float64[]

    for i in -25:25
        push!(ox, Float64(i))
        push!(oy, 15.0)
    end
    for i in -25:-4
        push!(ox, Float64(i))
        push!(oy, 4.0)
    end
    for i in -15:4
        push!(ox, -4.0)
        push!(oy, Float64(i))
    end
    for i in -15:4
        push!(ox, 4.0)
        push!(oy, Float64(i))
    end
    for i in 4:25
        push!(ox, Float64(i))
        push!(oy, 4.0)
    end
    for i in -4:4
        push!(ox, Float64(i))
        push!(oy, -15.0)
    end

    oox = ox[:]
    ooy = oy[:]

    sx = -10.0  # [m]
    sy = 6.0  # [m]
    syaw0 = deg2rad(00.0)
    syaw1 = deg2rad(00.0)

    @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    # sx = 14.0  # [m]
    # sy = 10.0  # [m]
    # syaw0 = deg2rad(00.0)
    # syaw1 = deg2rad(00.0)

    # @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    # sx = -14.0  # [m]
    # sy = 12.0  # [m]
    # syaw0 = deg2rad(00.0)
    # syaw1 = deg2rad(00.0)
    # @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    # sx = -20.0  # [m]
    # sy = 6.0  # [m]
    # syaw0 = deg2rad(00.0)
    # syaw1 = deg2rad(00.0)
    # @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    # sx = -14.0  # [m]
    # sy = 12.0  # [m]
    # syaw0 = deg2rad(00.0)
    # syaw1 = deg2rad(00.0)
    # path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    # sx = -20.0  # [m]
    # sy = 6.0  # [m]
    # syaw0 = deg2rad(180.0)
    # syaw1 = deg2rad(180.0)
    # @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    # sx = -20.0  # [m]
    # sy = 12.0  # [m]
    # syaw0 = deg2rad(180.0)
    # syaw1 = deg2rad(180.0)

    # @time path = calc_hybrid_astar_path(sx, sy, syaw0, syaw1, gx, gy, gyaw0, gyaw1, ox, oy, XY_GRID_RESOLUTION, YAW_GRID_RESOLUTION)

    # # Base.Test.@test length(path.x)>=1

    println("Test Done !!!")
end


if length(PROGRAM_FILE)!=0 &&
    contains(@__FILE__, PROGRAM_FILE)

    # test()
    main()
end


end #module

