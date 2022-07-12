import sys
import pymeshlab as ml
ms = ml.MeshSet()
ms.load_new_mesh(sys.argv[1])
ms.load_filter_script('meshlab_reset_origin.mlx')
ms.apply_filter_script()
ms.save_current_mesh(sys.argv[2])