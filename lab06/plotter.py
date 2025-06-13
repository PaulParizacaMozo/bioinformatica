import sys
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import matplotlib.cm as cm

def plot_dendrogram(linkage_file, labels_str, output_file, method_name):
    try:
        # carga la matriz
        Z = np.loadtxt(linkage_file)
        
        # carga etiquetas
        labels = labels_str.split(',')
        
        if Z.shape[0] + 1 != len(labels):
            print(f"Error de dimensiones: {Z.shape[0]+1} hojas vs {len(labels)} etiquetas.")
            sys.exit(1)

        # preparar la figura
        fig, ax = plt.subplots(figsize=(12, 8))
        plt.title(f'Dendrograma - Método: Distancia {method_name.capitalize()}', fontsize=16)
        plt.xlabel('Secuencias / Clusters', fontsize=12)
        plt.ylabel('Distancia', fontsize=12)
        
        # grafica el dendrograma
        ddata = dendrogram(
            Z,
            ax=ax,
            labels=labels,
            leaf_rotation=90.,
            leaf_font_size=10,
            link_color_func=lambda k: 'lightgrey'
        )

        # colores
        colors = cm.viridis(np.linspace(0, 1, len(ddata['icoord'])))
        
        for i, (icoords, dcoords) in enumerate(zip(ddata['icoord'], ddata['dcoord'])):
            color = colors[i]
            x1, y1 = icoords[0], dcoords[0]
            x2, y2 = icoords[1], dcoords[1]
            x3, y3 = icoords[2], dcoords[2]
            x4, y4 = icoords[3], dcoords[3]
            
            # Re-dibujar los segmentos de la unión con su color único
            ax.plot([x1, x2], [y1, y2], color=color, linewidth=1.5)
            ax.plot([x2, x3], [y2, y3], color=color, linewidth=1.5)
            ax.plot([x3, x4], [y3, y4], color=color, linewidth=1.5)

            x_center = 0.5 * (x2 + x3)
            y_center = y2
            ax.plot(x_center, y_center, 'o', color=color, markeredgewidth=0)

            ax.annotate(f"{y_center:.2f}", (x_center, y_center), xytext=(2, 2),
                        textcoords='offset points',
                        va='bottom', ha='left',
                        fontsize=8,
                        color='black')

        plt.grid(axis='y', linestyle='--')
        plt.tight_layout()
        plt.savefig(output_file)
        plt.close(fig)

    except FileNotFoundError:
        print(f"Error: No se encontró el archivo de enlace '{linkage_file}'.")
        sys.exit(1)
    except Exception as e:
        print(f"Ocurrió un error en el script de Python: {e}")
        sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python plotter.py <archivo_enlace> <etiquetas_csv> <archivo_salida_png> <nombre_metodo>")
        sys.exit(1)
    
    linkage_file = sys.argv[1]
    labels_str = sys.argv[2]
    output_file = sys.argv[3]
    method_name = sys.argv[4]
    
    plot_dendrogram(linkage_file, labels_str, output_file, method_name)
