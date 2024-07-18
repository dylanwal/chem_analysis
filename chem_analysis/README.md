
# High level structure of analysis


```mermaid
flowchart LR
    Signal --multiple--> Signal2D
    Signal2D --extract--> Signal
    Signal2D --multiple--> Signal3D
    Signal3D --extract--> Signal2D
```
* `unify_method` is used to convert mutiple signals into a 1D higher dimension. If the singals are all similar than the `unify_method` may only combine the data. However, if signals are all not uniform (same exactly the axis, data point for data point) then the `unify_method` may need to expand or cut or linearly interpolate to make a uniform axis. 


```mermaid
flowchart TB
    subgraph base
    Signal 
    Signal2D
    Signal3D
    end
 
    subgraph NMR
    Signal --> NMRSignal
    Signal2D --> NMRSignal2D
    end

    subgraph IR
    Signal --> IRSignal
    Signal2D --> IRSignal2D
    end
    
    subgraph MS
    Signal  --> MSSignal
    Signal2D --> MSSignal2D
    Signal3D --> MSSignal3D
    end
        
    subgraph SEC
    Signal --> SECSignal
    Signal2D --> SECSignal2D
    SECSignal --> SECCalibration 
    end
    
    subgraph gc_lc
        Signal --> GCSignal
        Signal --> GCMSSignal
        MSSignal2D -.-> GCMSSignal
        Signal2D --> GCMSSignal2D
        MSSignal3D -.-> GCMSSignal2D
    end

```

```mermaid
flowchart LR
    subgraph processing
        smoothing --> baseline_corrections
    end

    subgraph analysis
        peak_picking --> boundary_detection --> integration
    end

    Signal --> processing --> analysis
```

## Important Notes
- When indexing; the first index referses to the outer value
- Data is always stored in sorted order (small to large)
  - Adjust order in plotting, if data is typically presented in reverse.