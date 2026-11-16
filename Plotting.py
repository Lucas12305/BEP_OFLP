def plot_solution(OFP, instance : dict, save_path : str = None) -> None:
    '''
    This makes a plot of the solution that SCIP has found for the obnoxious facility location problem.
    Around every community it draws a grey circle with the radius of the objective, the facilities are placed outside of the grey area.
    :param OFP: The model that has been optimized.
    :param instance: The instance that has been solved, as given by read_instance.
    :param save_path: If this is given the plot is saved to this file instead of shown.
    '''
    solution = get_solution(OFP, instance)
    if solution is None:
        print("There is no solution to plot.")
        return None
    facilities, L = solution

    x_lb, x_ub = instance['x_bounds']
    y_lb, y_ub = instance['y_bounds']
    communities = instance['communities']

    # Create figure and axes
    fig, ax = plt.subplots()

    # Create a geometric box (Rectangle)
    box = patches.Rectangle((x_lb, y_lb), x_ub - x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)

    # Plot the communities with a grey circle around them, the facilities have to be outside of the grey area
    # L is a squared distance, so the radius is the square root of L
    # The fill is see-through, so overlapping circles get darker and every circle stays visible
    for c in communities:
        circle = patches.Circle((c[0], c[1]), np.sqrt(L), linewidth=0.8, edgecolor='gray', facecolor=(0.5, 0.5, 0.5, 0.15))
        ax.add_patch(circle)
        ax.scatter(c[0], c[1], s=10, c='black', marker='o', zorder=20)

    # Plot the facilities
    for i, f in enumerate(facilities):
        ax.scatter(f[0], f[1], s=70, c='red', zorder=21)
        ax.annotate(f"F{i}", (f[0], f[1]), textcoords="offset points", xytext=(6, 6), zorder=22)

    # Set limits and equal aspect ratio
    ax.set_xlim(x_lb - abs(x_lb - x_ub) / 2, x_ub + abs(x_lb - x_ub) / 2)
    ax.set_ylim(y_lb - abs(y_lb - y_ub) / 2, y_ub + abs(y_lb - y_ub) / 2)
    ax.set_aspect('equal')
    ax.set_title(f"L = {L:.4f}, radius = {np.sqrt(L):.4f}")

    # Show the plot, or save it when a file is given
    plt.grid(True)
    if save_path is None:
        plt.show()
    else:
        fig.savefig(save_path, dpi=150)
        plt.close(fig)
        print(f"The plot has been saved to {save_path}")
    return None





def plot_box_LP(x_lb,x_ub,y_lb,y_ub, circles, intersections, planes : list, trc :bool = False, rc : list=[]) -> None:
    '''
    This makes a plot for the box that gets cut using the LP.
    '''
    # Create figure and axes
    fig, ax = plt.subplots()

    # Create a geometric box (Rectangle)
    # Arguments: (lower-left x, y), width, height
    box = patches.Rectangle((x_lb, y_lb), x_ub-x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)
        
    # Create the circles
    # Arguments: (center x, y), radius
    for c in circles:
        circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='black', fill = False)
        ax.add_patch(circle)
    # Add shapes to the axes

    for p in intersections:
        ax.scatter(p[0], p[1], s=20, c='red', marker='o')
    if trc == True:
        for c in rc:
            circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='lightgray', fill = False)
            ax.add_patch(circle)
            ax.scatter(c[0], c[1], s=1, c='blue', marker='o', label='Center Point')
        
    for plane in planes:
        p0, p1, p2 = plane
    
        if abs(p1) > 1e-10:
            x_vals = np.linspace(x_lb, x_ub, 200)
            y_line = -(p0 * x_vals + p2) / p1
    
            # Keep everything inside the box
            y_line = np.clip(y_line, y_lb, y_ub)
    
            if p1 < 0:
                # Shade below line
                ax.fill_between(
                    x_vals, y_lb, y_line,
                    alpha=1, color = 'blue'
                )
            else:
                # Shade above line
                ax.fill_between(
                    x_vals, y_line, y_ub,
                    alpha=1, color = 'blue'
                )
    
            ax.plot(x_vals, y_line, color = 'blue')
    
        else:
            # Vertical line
            if abs(p0) > 1e-10:
                x_line = -p2 / p0
    
                if x_lb <= x_line <= x_ub:
                    ax.vlines(x_line, y_lb, y_ub, color='blue')
    
                    if p0 < 0:
                        # x <= x_line
                        ax.axvspan(x_lb, x_line, alpha=1, color = 'blue')
                    else:
                        # x >= x_line
                        ax.axvspan(x_line, x_ub, alpha=1, color = 'blue')
        
    # Set limits and equal aspect ratio
    ax.set_xlim(x_lb - abs(x_lb-x_ub)/2, x_ub + abs(x_lb-x_ub)/2)
    ax.set_ylim(y_lb - abs(y_lb-y_ub)/2, y_ub + abs(y_lb-y_ub)/2)
    ax.set_aspect('equal')
        
    # Show the plot
    plt.grid(True)
    plt.show()
    return None





def plot_box_GS(x_lb,x_ub,y_lb,y_ub, circles, intersections, sorted_points : list, trc :bool = False, rc : list=[]) -> None:
    '''
    This makes a plot for the box that gets cut using the LP.
    '''
    # Create figure and axes
    fig, ax = plt.subplots()

    # Create a geometric box (Rectangle)
    # Arguments: (lower-left x, y), width, height
    box = patches.Rectangle((x_lb, y_lb), x_ub-x_lb, y_ub - y_lb, linewidth=3, edgecolor='black', facecolor='none')
    ax.add_patch(box)
        
    # Create the circles
    # Arguments: (center x, y), radius
    for c in circles:
        circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='black', fill = False)
        ax.add_patch(circle)
    
    # Add shapes to the axes
    for p in intersections:
        #print(f"We have a point {p}")
        ax.scatter(p[0], p[1], s=70, c='green', zorder=21)
    
    if trc == True:
        for c in rc:
            circle = patches.Circle((c[0], c[1]), c[2], linewidth=1, edgecolor='lightgray', fill = False)
            ax.add_patch(circle)
            ax.scatter(c[0], c[1], s=1, c='blue', marker='o', label='Center Point')
        
    for i in range(len(sorted_points)):
        p1 = sorted_points[i]
        p2 = sorted_points[(i + 1) % len(sorted_points)]
        ax.plot(
            [p1[0], p2[0]],
            [p1[1], p2[1]],
            color='blue'
        )
  
    # Set limits and equal aspect ratio
    ax.set_xlim(x_lb - abs(x_lb-x_ub)/2, x_ub + abs(x_lb-x_ub)/2)
    ax.set_ylim(y_lb - abs(y_lb-y_ub)/2, y_ub + abs(y_lb-y_ub)/2)
    ax.set_aspect('equal')
        
    # Show the plot
    plt.grid(True)
    plt.show()
    return None