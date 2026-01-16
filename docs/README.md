# Documentación técnica

## Estrategia para GNINA

- **Disponibilidad en conda-forge/bioconda:** no existe el paquete `gnina` en esos canales (verificado vía Anaconda API).
- **Estrategia aplicada:** el entorno conda usa `smina` para tareas compatibles, y el rescoring con GNINA se ejecuta mediante el contenedor `gnina_container_wrapper.sh`. Esta separación evita depender de un paquete inexistente en los canales principales.

### Archivos relacionados

- `workflow/envs/gnina.yaml`: define el entorno conda con `smina`.
- `scripts/gnina_container_wrapper.sh`: wrapper para ejecutar GNINA en contenedor.
